#include "JetMatchingSubjets.h"
#include "ZsjTools.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <g4histos/G4VtxNtuple.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>
#include <TVector2.h>  // for TVector2::Phi_mpi_pi
#include <jetbase/JetMap.h>

#include <centrality/CentralityInfo.h>
#include <jetbackground/TowerBackground.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <TH1F.h>
#include <TH1D.h>
#include <TTree.h>
#include <RooUnfoldResponse.h>

#include <TParameter.h>
#include <TNamed.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace {
  template<class T> T* book1D(const char* name, const char* title, const std::vector<double>& edges) {
    return new T(name, title, static_cast<int>(edges.size()-1), edges.data());
  }
}

JetMatchingSubjets::JetMatchingSubjets(const std::string& recojetname,
                                       const std::string& truthjetname,
                                       const std::string& outputfilename)
: SubsysReco(std::string("JetMatchingSubjets_")+recojetname+"_"+truthjetname)
, m_recoJetName(recojetname)
, m_truthJetName(truthjetname)
, m_outputFileName(outputfilename)
{}

JetMatchingSubjets::~JetMatchingSubjets() {}

int JetMatchingSubjets::Init(PHCompositeNode*)
{
  if (m_ptEdges.empty()) m_ptEdges = {5,10,15,20,25,30,35,40,45,50,55};
  m_nPtBins = (int)m_ptEdges.size() - 1;
  
  if (m_zsjEdges.empty()) { const int n=10; for (int i=0;i<=n;++i) m_zsjEdges.push_back(0.5*i/n); }
  m_nZsjBins = (int)m_zsjEdges.size() - 1;
  
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  
  //TTree building
  m_T = new TTree("T", "Jet matching for RooUnfold (pt and pt×z_sj)");
  m_T->Branch("event", &m_event, "event/I");
  m_T->Branch("cent",  &m_centrality, "cent/F");
  m_T->Branch("b",     &m_b, "b/F");
  //m_T->Branch("vtx_z", &m_vtx_z);
  m_T->Branch("vtx_z", &m_vtx_z, "vtx_z/F");
  m_T->Branch("reco_pt",  &v_reco_pt);
  m_T->Branch("reco_eta", &v_reco_eta);
  m_T->Branch("reco_phi", &v_reco_phi);
  m_T->Branch("reco_zsj", &v_reco_zsj);
  m_T->Branch("reco_thetasj", &v_reco_thetasj);
  
  m_T->Branch("truth_pt",  &v_truth_pt);
  m_T->Branch("truth_eta", &v_truth_eta);
  m_T->Branch("truth_phi", &v_truth_phi);
  m_T->Branch("truth_zsj", &v_truth_zsj);
  m_T->Branch("truth_thetasj", &v_truth_thetasj);
  
  // matched pairs
  m_T->Branch("match_reco_pt",  &v_match_reco_pt);
  m_T->Branch("match_reco_eta", &v_match_reco_eta);
  m_T->Branch("match_reco_phi", &v_match_reco_phi);
  m_T->Branch("match_reco_zsj", &v_match_reco_zsj);
  m_T->Branch("match_reco_thetasj", &v_match_reco_thetasj);
  
  m_T->Branch("match_truth_pt",  &v_match_truth_pt);
  m_T->Branch("match_truth_eta", &v_match_truth_eta);
  m_T->Branch("match_truth_phi", &v_match_truth_phi);
  m_T->Branch("match_truth_zsj", &v_match_truth_zsj);
  m_T->Branch("match_truth_thetasj", &v_match_truth_thetasj);
  m_T->Branch("match_dR",        &v_match_dR);

  // fakes & misses
  m_T->Branch("fake_reco_pt",  &v_fake_reco_pt);
  m_T->Branch("fake_reco_eta", &v_fake_reco_eta);
  m_T->Branch("fake_reco_phi", &v_fake_reco_phi);
  m_T->Branch("fake_reco_zsj", &v_fake_reco_zsj);
  m_T->Branch("fake_reco_thetasj", &v_fake_reco_thetasj);
  
  m_T->Branch("fake_truth_pt",  &v_fake_truth_pt);
  m_T->Branch("fake_truth_eta", &v_fake_truth_eta);
  m_T->Branch("fake_truth_phi", &v_fake_truth_phi);
  m_T->Branch("fake_truth_zsj", &v_fake_truth_zsj);
  m_T->Branch("fake_truth_thetasj", &v_fake_truth_thetasj);
  
  return Fun4AllReturnCodes::EVENT_OK;
}
static inline double lookupZsjForReco(
				      const Jet* r,
				      const std::unordered_map<const Jet*, double>& zReco,
				      double dr_tol = 0.2, double pt_tol = 0.5)
{
  auto it = zReco.find(r);
  if (it != zReco.end()) return it->second;
  
  // Pointer aliasing fallback: match by ΔR and pT
  const double eta = r->get_eta();
  const double phi = r->get_phi();
  double best_dr2 = 1e9; //start with something impossibly large such that the iterations can become smaller and smaller
  double zr = std::numeric_limits<double>::quiet_NaN();
  
  for (const auto& kv : zReco) {
    const Jet* rc = kv.first;
    const double dphi = TVector2::Phi_mpi_pi(phi - rc->get_phi());
    const double deta = eta - rc->get_eta();
    const double dr2  = dphi*dphi + deta*deta;
    if (dr2 <= dr_tol*dr_tol && std::abs(r->get_pt() - rc->get_pt()) <= pt_tol) {
      if (dr2 < best_dr2) { best_dr2 = dr2; zr = kv.second; }
    }
  }
  return zr; // may remain NaN if no match
}

static Jet* matchRaw(const Jet* jcal, JetContainer* raw, const double dRmax=5e-2)
{
  if (!jcal || !raw) return nullptr;
  Jet* best = nullptr;
  double bestDR2 = dRmax*dRmax;
  const double eta = jcal->get_eta();
  const double phi = jcal->get_phi();
  
  for (auto it = raw->begin(); it != raw->end(); ++it) {
    Jet* jr = *it;
    const double dphi = TVector2::Phi_mpi_pi(phi - jr->get_phi());
    const double deta = eta - jr->get_eta();
    const double dr2  = dphi*dphi + deta*deta;
    if (dr2 < bestDR2) { bestDR2 = dr2; best = jr; }
  }
  return best;
}

int JetMatchingSubjets::process_event(PHCompositeNode* topNode)
{
  ++m_event;
  int nRecoNoRaw = 0;
  
  // clear vectors
  v_reco_pt.clear();  v_reco_eta.clear();  v_reco_phi.clear();  v_reco_zsj.clear(); v_reco_thetasj.clear();
  v_truth_pt.clear(); v_truth_eta.clear(); v_truth_phi.clear(); v_truth_zsj.clear(); v_truth_thetasj.clear();
  v_match_reco_pt.clear(); v_match_reco_eta.clear(); v_match_reco_phi.clear(); v_match_reco_zsj.clear(); v_match_reco_thetasj.clear();
  v_match_truth_pt.clear(); v_match_truth_eta.clear(); v_match_truth_phi.clear(); v_match_truth_zsj.clear(); v_match_truth_thetasj.clear();
  v_match_dR.clear();
  v_fake_reco_pt.clear(); v_fake_reco_eta.clear(); v_fake_reco_phi.clear(); v_fake_reco_zsj.clear(); v_fake_reco_thetasj.clear();
  v_fake_truth_pt.clear(); v_fake_truth_eta.clear(); v_fake_truth_phi.clear(); v_fake_truth_zsj.clear(); v_fake_truth_thetasj.clear();

  recoToTruth.clear();
  truthToReco.clear();
  
  auto* jetsReco  = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  auto* jetsTruth = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  auto* cent      = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  auto* em   = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  auto* ih   = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
  auto* oh   = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
  auto* geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  auto* geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto* geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  auto* bg    = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub1");

  if (!jetsReco || !jetsTruth || !cent || !em || !ih || !oh || !geomEM || !geomIH || !geomOH || !bg) {
    std::cerr << "[JetMatchingSubjets] Missing node(s)." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_centrality = cent->get_centile(CentralityInfo::PROP::bimp);
  m_b          = cent->get_quantity(CentralityInfo::PROP::bimp);

  // Truth vertex (simulation)
  auto* truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (truthInfo) {
    const PHG4VtxPoint* vtx = truthInfo->GetPrimaryVtx(truthInfo->GetPrimaryVertexIndex());
    m_vtx_z = vtx ? static_cast<float>(vtx->get_z()) : std::numeric_limits<float>::quiet_NaN();
  } else {
    m_vtx_z = std::numeric_limits<float>::quiet_NaN();
  }
  
  // Enforce ±30 cm *before* any jet work. This keeps only in-window events in the TTree.
  if (!std::isfinite(m_vtx_z) || m_vtx_z < m_zmin || m_vtx_z >= m_zmax) {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  ++m_N_in_zwin;
  auto* jetsRecoCal = findNode::getClass<JetContainer>(topNode, m_recoJetName.c_str());
  JetContainer* jetsRecoRaw = jetsRecoCal;
  if (!m_recoConstituentJetNode.empty()) {
    jetsRecoRaw = findNode::getClass<JetContainer>(topNode, m_recoConstituentJetNode.c_str());
  }
  
  if (!jetsRecoCal) {
    std::cerr << "[JetMatchingSubjets] Missing calibrated jet node: " << m_recoJetName << "\n";
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!jetsRecoRaw) {
    std::cerr << "[JetMatchingSubjets] Missing raw constituent jet node: " << m_recoConstituentJetNode << "\n";
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  //  std::unordered_map<const Jet*, double> zReco;
  std::unordered_map<const Jet*, double> zReco;
  std::unordered_map<const Jet*, double> zTruth;
  std::unordered_map<const Jet*, double> thetaReco;
  std::unordered_map<const Jet*, double> thetaTruth;
  
  // Matching
  //MatchJets1to1(jetsReco, jetsTruth);
  MatchJets1to1(jetsRecoCal, jetsTruth);
  // try to compute z_sj/theta_sj from RAW twin

  // --- Reco loop: compute once on RAW constituents, store, and fill "reco_*" ---
  for (auto it = jetsRecoCal->begin(); it != jetsRecoCal->end(); ++it) {
    Jet* rcal = *it;
    if (rcal->get_pt() < m_recoPtMin) continue;
    if (std::abs(rcal->get_eta()) > m_etaRange.second) continue;
    
    // always push pt/eta/phi once per jet
    v_reco_pt .push_back(rcal->get_pt());
    v_reco_eta.push_back(rcal->get_eta());
    v_reco_phi.push_back(rcal->get_phi());
    double z = std::numeric_limits<float>::quiet_NaN();
    double th = std::numeric_limits<float>::quiet_NaN();
    
    Jet* rraw = (jetsRecoRaw == jetsRecoCal) ? rcal : matchRaw(rcal, jetsRecoRaw, /*dRmax=*/0.05);
    if (!rraw) {
      ++nRecoNoRaw;                          // count missing raw twin
    } else {
      double ztmp, thtmp;
      bool ok = ComputeZsjForJet(/*jet=*/rraw, em, ih, oh, geomEM, geomIH, geomOH, bg,
				 /*v2=*/0, /*psi=*/0, /*doUnsub=*/false,
				 /*R_jet=*/0.4, /*R_sub=*/0.2, /*ptCutSubjet=*/0.1,
				 ztmp, thtmp);
      if (ok && std::isfinite(ztmp) && std::isfinite(thtmp)) {
	z = ztmp;
	th = thtmp;
      }
    }
    
    // always cache something (NaN is fine if it failed)
    zReco[rcal]     = z;
    thetaReco[rcal] = th;
    
    // and still push to your per-jet vectors:
    v_reco_zsj.push_back(static_cast<float>(z));
    v_reco_thetasj.push_back(static_cast<float>(th));
  }
  std::cout << "[DBG] after reco cache: zReco.size=" << zReco.size()
	    << " thetaReco.size=" << thetaReco.size() << '\n';

  // Truth jets
  for (auto j: *jetsTruth) {
    if (j->get_pt() < m_truthPtMin) continue;
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;
    
    v_truth_pt.push_back(j->get_pt());
    v_truth_eta.push_back(j->get_eta());
    v_truth_phi.push_back(j->get_phi());

    double z, th;
    if (ComputeZsjForTruthJet(j, truthInfo, /*R_jet=*/0.4, /*R_sub=*/0.2, /*ptCutSubjet=*/0.1, z, th)) {
      zTruth[j] = z;
      thetaTruth[j] = th;
      
      v_truth_zsj.push_back(z);
      std::cout << "truth z_sj = " << z << std::endl;   // (optional debug)
      v_truth_thetasj.push_back(th);
      std::cout << "truth theta_sj = " << th << std::endl;   // (optional debug)
    } else {
      v_truth_zsj.push_back(NAN);
      v_truth_thetasj.push_back(NAN);
    } 
  }

  // Misses (truth-only)
  for (auto t: *jetsTruth) {
    if (t->get_pt() < m_truthPtMin) continue;
    if (std::abs(t->get_eta()) > m_etaRange.second) continue;
    if (truthToReco.find(t) != truthToReco.end()) continue;

    double zt = zTruth.count(t) ? zTruth[t] : NAN;
    double tht= thetaTruth.count(t) ? thetaTruth[t] : NAN;
    
    v_fake_truth_pt.push_back(t->get_pt());
    v_fake_truth_eta.push_back(t->get_eta());
    v_fake_truth_phi.push_back(t->get_phi());
    v_fake_truth_zsj.push_back(zt);
    v_fake_truth_thetasj.push_back(tht);
  }

  // --- Matched ---
  for (auto& kv: recoToTruth) {
    Jet* r = kv.first;
    Jet* t = kv.second;
    
    const float dr = DeltaR(r,t);
    v_match_reco_pt .push_back(r->get_pt());
    v_match_reco_eta.push_back(r->get_eta());  // FIXED: eta
    v_match_reco_phi.push_back(r->get_phi());  // FIXED: phi
    v_match_dR.push_back(dr);
    
    // truth side unchanged
    v_match_truth_pt .push_back(t->get_pt());
    v_match_truth_eta.push_back(t->get_eta());
    v_match_truth_phi.push_back(t->get_phi());
    double zt = zTruth.count(t) ? zTruth[t] : std::numeric_limits<double>::quiet_NaN();
    v_match_truth_zsj.push_back(zt);
    double tht = thetaTruth.count(t) ? thetaTruth[t] : std::numeric_limits<double>::quiet_NaN();
    v_match_truth_thetasj.push_back(tht);
    
    // inside matched loop, before lookup
    auto it = zReco.find(r);
    std::cout << "r@" << r << " in-cache? " << (it!=zReco.end()) << "   cache size=" << zReco.size() << '\n';
    
    double zr = lookupZsjForReco(r, zReco, 0.01, 0.1);  // no third argument now
    v_match_reco_zsj.push_back(std::isfinite(zr) ? zr : std::numeric_limits<float>::quiet_NaN());
    double thr = lookupZsjForReco(r, thetaReco, 0.01, 0.1);
    v_match_reco_thetasj.push_back(std::isfinite(thr) ? thr : std::numeric_limits<float>::quiet_NaN());
    // optional debug
    std::cout << "match truth z_sj = " << zt << "\nmatch reco z_sj = " << zr << std::endl;
    std::cout << "match truth theta_sj = " << tht << "\nmatch reco theta_sj = " << thr << std::endl;
    
  }
  
  // --- Fakes (reco-only): iterate the SAME calibrated container used to build zReco ---
  std::unordered_set<const Jet*> matchedReco;
  matchedReco.reserve(recoToTruth.size());
  int nLookupNaN = 0;
  for (auto& kv : recoToTruth) matchedReco.insert(kv.first);
  
  for (auto r: *jetsRecoCal) {   // <- use jetsRecoCal, not jetsReco
    if (r->get_pt() < m_recoPtMin) continue;
    if (std::abs(r->get_eta()) > m_etaRange.second) continue;
    if (matchedReco.count(r)) continue;  // robust skip
    
    v_fake_reco_pt .push_back(r->get_pt());
    v_fake_reco_eta.push_back(r->get_eta());
    v_fake_reco_phi.push_back(r->get_phi());
    
    //double zr = lookupZsjForReco(r, zReco, DeltaR);
    double zr = lookupZsjForReco(r, zReco, 0.01, 0.1);  // no third argument now
    v_fake_reco_zsj.push_back(std::isfinite(zr) ? zr : std::numeric_limits<float>::quiet_NaN());
    double thr = lookupZsjForReco(r, thetaReco, 0.01, 0.1);
    v_fake_reco_thetasj.push_back(std::isfinite(thr) ? thr : std::numeric_limits<float>::quiet_NaN());
    if (std::isnan(zr) || std::isnan(thr)) {
      ++nLookupNaN;
    }
  }
  std::cout << "[Evt " << m_event << "] "
	    << "RecoNoRaw=" << nRecoNoRaw	    
	    << " LookupNaN="<< nLookupNaN
          << std::endl;
  
  m_T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}
int JetMatchingSubjets::End(PHCompositeNode*)
{
  PHTFileServer::get().cd(m_outputFileName);
  if (m_T) { m_T->Write(); m_T = nullptr; }
  // Store N(z-win) and sigma for this file

  TParameter<double>    pS("sigma_pb",  m_sigma_pb);
  TNamed                pName("sample_name", m_sample_name.c_str());  
  TParameter<long long>("N_in_zwin", m_N_in_zwin).Write(); // number of entries written (|z|<=30)
  pS.Write();
  pName.Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::Reset(PHCompositeNode*) { return Fun4AllReturnCodes::EVENT_OK; }
void JetMatchingSubjets::Print(const std::string& what) const { }
  
// ---- helper definitions ----
static inline float wrapPhi(float dphi){
  while (dphi >  M_PI) dphi -= 2*M_PI;
  while (dphi < -M_PI) dphi += 2*M_PI;
  return dphi; 
}

float JetMatchingSubjets::DeltaR(Jet* a, Jet* b)
{
  const float dEta = a->get_eta() - b->get_eta(); 
  float dPhi = wrapPhi(a->get_phi() - b->get_phi());
  return std::sqrt(dEta*dEta + dPhi*dPhi);
}

void JetMatchingSubjets::MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets)
{
  struct Pair { float dr; Jet* r; Jet* t; };  
  std::vector<Pair> cand;
  for (auto r: *recoJets) {
    if (r->get_pt() < m_recoPtMin) continue;
    if (std::abs(r->get_eta()) > m_etaRange.second) continue;
    for (auto t: *truthJets) {
      if (t->get_pt() < m_truthPtMin) continue;
      if (std::abs(t->get_eta()) > m_etaRange.second) continue;
      float dr = DeltaR(r,t);
      if (dr < m_matchDRMax) cand.push_back({dr, r, t});
    }
  }
  std::sort(cand.begin(), cand.end(), [](const Pair& a, const Pair& b){ return a.dr < b.dr; });
  std::set<Jet*> usedR, usedT;
  recoToTruth.clear();
  truthToReco.clear();
  for (auto& p: cand) {
    if (usedR.count(p.r) || usedT.count(p.t)) continue;
    usedR.insert(p.r);
    usedT.insert(p.t);
    recoToTruth[p.r] = p.t;
    truthToReco[p.t] = p.r;
     
  }
}

