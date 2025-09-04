#ifndef JETMATCHING_H
#define JETMATCHING_H

#include <fun4all/SubsysReco.h>
#include <limits>
#include <string>
#include <vector>
#include <utility>
#include <set>
#include <map>
#include <functional>

class PHCompositeNode;
class TTree;
class TH1F;
class TH1D;
class TH2D;
class RooUnfoldResponse;
class Jet;
class JetContainer;
class CentralityInfo;

class JetMatching : public SubsysReco
{
public:
  JetMatchingSubjets(const std::string& recojetname,
                     const std::string& truthjetname,
                     const std::string& outputfilename);
  ~JetMatchingSubjets() override;

  // Config
  void setEtaRange(double low, double high) { m_etaRange = {low, high}; }
  void setRecoPtMin(double x) { m_recoPtMin = x; }
  void setTruthPtMin(double x) { m_truthPtMin = x; }
  void setMatchDRMax(double x) { m_matchDRMax = x; }
  void setPtBinning(const std::vector<double>& edges) { m_ptEdges = edges; }
  void setZWindow(float zmin, float zmax) { m_zmin = zmin; m_zmax = zmax; }
  void setRecoConstituentJetNode(const std::string& n) { m_recoConstituentJetNode = n; }
  
  // Fun4All
  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  void Print(const std::string& what) const override;

private:
  // Helpers
  void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets);
  static float DeltaR(Jet* a, Jet* b);
  
private:
  // Names / IO
  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;
  std::string m_recoConstituentJetNode;
 
  // Kinematic selections
  std::pair<double,double> m_etaRange { -0.7, 0.7 };
  double m_recoPtMin  = 5.0;
  double m_truthPtMin = 10.0;
  double m_matchDRMax = 0.2;

  // Binning
  std::vector<double> m_ptEdges;   // default set in Init()
  int m_nPtBins  = 0;

  // Matching maps (1↔1 greedy)
  std::map<Jet*, Jet*> recoToTruth;
  std::map<Jet*, Jet*> truthToReco;

  // Tree
  TTree* m_T = nullptr;
  int   m_event = -1;
  float m_centrality = 0;
  float m_b = 0;
  // --- z-vertex and metadata for per-file bookkeeping ---
  float      m_vtx_z    = std::numeric_limits<float>::quiet_NaN();
  long long  m_N_in_zwin = 0;
  
  // configurable z-window (used only if you want to count inside the module)
  float      m_zmin = -30.f;
  float      m_zmax =  30.f;
  
  // optional: store per-file sigma and sample name (purely for provenance)
  double     m_sigma_pb = 0.0;
  std::string m_sample_name = "unknown";
  
  void setSigmaPb(double s) { m_sigma_pb = s; }
  void setSampleName(const std::string& n) { m_sample_name = n; }
  
  // Per-event containers
  std::vector<float> v_reco_pt,  v_reco_eta,  v_reco_phi;
  std::vector<float> v_truth_pt, v_truth_eta, v_truth_phi;

  // Matched pairs (parallel arrays, same length)
  std::vector<float> v_match_reco_pt,  v_match_reco_eta,  v_match_reco_phi;
  std::vector<float> v_match_truth_pt, v_match_truth_eta, v_match_truth_phi;
  std::vector<float> v_match_dR;

  // Fakes / Misses
  std::vector<float> v_fake_reco_pt,  v_fake_reco_eta,  v_fake_reco_phi;   // detector-only (no match)
  std::vector<float> v_fake_truth_pt, v_fake_truth_eta, v_fake_truth_phi;  // truth-only (missed)
};

#endif
