
#include "ZsjTools.h"

#include <jetbase/Jet.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerDefs.h>
#include <jetbackground/TowerBackground.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <cmath>

using namespace fastjet;
using fastjet::contrib::SoftDrop;
std::vector<PseudoJet>
BuildPseudoJets(Jet* jet,
                TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
                TowerBackground* bg, float v2, float psi2, bool doUnsub)
{
  std::vector<PseudoJet> particles;
  if (!jet) return particles;
  static int dbg_prints = 0;
  int nTried = 0, nTower = 0, nGeo = 0, nEpos = 0;
  int nCEMC = 0, nIH = 0, nOH = 0;
  
  for (const auto& comp : jet->get_comp_vec())
  {
    const unsigned int ch = comp.second;

    TowerInfo* tower = nullptr;                 // <- define ONCE
    const RawTowerGeom* geo = nullptr;          // <- define ONCE
    unsigned int calokey = 0;
    int ieta = 0, iphi = 0;
    float eta = 0.f, phi = 0.f, UE = 0.f;
    ++nTried;
    if (comp.first == 14 || comp.first == 29) {
      ++nCEMC;
      // -------- CEMC --------
      if (!em || !geomEM) continue;
      tower = em->get_tower_at_channel(ch);
      if (tower) {
        calokey = em->encode_key(ch);           // channel -> key
	++nTower;
      } else {
        // If your TowerInfoContainer lacks get_tower_at_key, remove this block.
        tower = em->get_tower_at_key(ch);       // try treating ch as a key
        if (!tower) continue;
        calokey = ch;                           // ch is already the key
      }
      ieta = em->getTowerEtaBin(calokey);
      iphi = em->getTowerPhiBin(calokey);
      const auto key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, ieta, iphi);
      geo = geomEM->get_tower_geometry(key);
      if (!geo) continue;
      ++nGeo;
      eta = geo->get_eta();
      phi = geo->get_phi();
      UE  = (bg ? bg->get_UE(0).at(ieta) : 0.f);

    } else if (comp.first == 15 || comp.first == 30) {
      ++nIH;
      // -------- HCALIN --------
      if (!ih || !geomIH) continue;

      tower = ih->get_tower_at_channel(ch);
      if (tower) {
        calokey = ih->encode_key(ch);
	++nTower;
      } else {
        // If not available in your build, remove this fallback.
        tower = ih->get_tower_at_key(ch);
        if (!tower) continue;
        calokey = ch;
      }

      ieta = ih->getTowerEtaBin(calokey);
      iphi = ih->getTowerPhiBin(calokey);
      const auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALIN, ieta, iphi);
      geo = geomIH->get_tower_geometry(key);
      if (!geo) continue;
      ++nGeo;
      eta = geo->get_eta();
      phi = geo->get_phi();
      UE  = (bg ? bg->get_UE(1).at(ieta) : 0.f);

    } else if (comp.first == 16 || comp.first == 31) {
      ++nOH;
      // -------- HCALOUT --------
      if (!oh || !geomOH) continue;

      tower = oh->get_tower_at_channel(ch);
      if (tower) {
        calokey = oh->encode_key(ch);
	++nTower;
      } else {
        // If not available in your build, remove this fallback.
        tower = oh->get_tower_at_key(ch);
        if (!tower) continue;
        calokey = ch;
      }

      ieta = oh->getTowerEtaBin(calokey);
      iphi = oh->getTowerPhiBin(calokey);
      const auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALOUT, ieta, iphi);
      geo = geomOH->get_tower_geometry(key);
      if (!geo) continue;
      ++nGeo;
      eta = geo->get_eta();
      phi = geo->get_phi();
      UE  = (bg ? bg->get_UE(2).at(ieta) : 0.f);

    } else {
      continue; // unknown comp type
    }

    // UE modulation AFTER we know phi
    const float UE_mod = UE * (1.f + 2.f * v2 * std::cos(2.f * (phi - psi2)));

    // tower was validated above; still, be defensive
    if (!tower) continue;
    float energy = tower->get_energy();
    if (doUnsub) energy -= UE_mod;

    // drop non-physical or over-subtracted towers
    if (!std::isfinite(energy) || energy <= 0.f) continue;

    const float coshEta = std::cosh(eta);
    if (!(coshEta > 0.f)) continue;

    const float pt = energy / coshEta;
    if (!std::isfinite(pt) || pt <= 0.f) continue;

    const float px = pt * std::cos(phi);
    const float py = pt * std::sin(phi);
    const float pz = pt * std::sinh(eta);

    particles.emplace_back(px, py, pz, energy);
  }
  if (particles.empty() && dbg_prints < 12) {
    std::cout << "[DBG] BuildPseudoJets reco: tried=" << nTried
	      << " tower=" << nTower << " geo=" << nGeo
	      << " E>0=" << nEpos
	      << " | CEMC=" << nCEMC << " IH=" << nIH << " OH=" << nOH << "\n";
    ++dbg_prints;
  }
  
  return particles;
}

bool ComputeZsjForJet(Jet* jet,
                      TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                      RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
                      TowerBackground* bg, float v2, float psi2,
                      bool doUnsub,
                      double R_jet, double R_sub, double ptCutSubjet,
                      double& zsj, double& theta)
{
  zsj = NAN; theta = NAN;
  auto particles = BuildPseudoJets(jet, em, ih, oh, geomEM, geomIH, geomOH, bg, v2, psi2, doUnsub);
  if (particles.empty()) return false;

  fastjet::JetDefinition jetDefAKT_R(fastjet::antikt_algorithm, R_jet);
  fastjet::ClusterSequence clustSeq(particles, jetDefAKT_R);
  auto jets = fastjet::sorted_by_pt(clustSeq.inclusive_jets());
  if (jets.empty()) return false;

  fastjet::PseudoJet leading = jets[0];
  // Subjet reclustering
  fastjet::JetDefinition jetDefSub(fastjet::antikt_algorithm, R_sub);
  fastjet::ClusterSequence subClust(leading.constituents(), jetDefSub);
  auto subjets = fastjet::sorted_by_pt(subClust.inclusive_jets());
  if (subjets.size() < 2) return false;

  auto sj1 = subjets[0];
  auto sj2 = subjets[1];
  if (sj1.perp() < ptCutSubjet || sj2.perp() < ptCutSubjet) return false;

  theta = sj1.delta_R(sj2);
  zsj   = sj2.perp() / (sj1.perp() + sj2.perp());
  if (!std::isfinite(theta) || !std::isfinite(zsj)) return false;
  return true;
}

bool ComputeZsjForTruthJet(Jet* truthJet,
                           PHG4TruthInfoContainer* truthInfo,
                           double R_jet, double R_sub, double ptCutSubjet,
                           double& zsj, double& theta)
{
  zsj = NAN; theta = NAN;
  if (!truthJet || !truthInfo) return false;

  std::vector<fastjet::PseudoJet> parts;
  parts.reserve(truthJet->size_comp());

  // Rebuild constituents from PHG4 particles referenced by this truth jet
  for (const auto& comp : truthJet->get_comp_vec()) {
    // comp.second is typically the track id for truth-particle constituents
    const int trackid = static_cast<int>(comp.second);
    const PHG4Particle* p = truthInfo->GetParticle(trackid);
    if (!p) continue;

    const double px = p->get_px();
    const double py = p->get_py();
    const double pz = p->get_pz();
    const double E  = p->get_e();
    // require positive energy and finite kinematics
    if (!(std::isfinite(px) && std::isfinite(py) && std::isfinite(pz) && std::isfinite(E))) continue;
    if (E <= 0.0) continue;

    parts.emplace_back(px, py, pz, E);
  }
  if (parts.size() < 2) return false;

  // cluster jet and subjets
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R_jet);
  fastjet::ClusterSequence cs(parts, jetDef);
  auto jets = fastjet::sorted_by_pt(cs.inclusive_jets());
  if (jets.empty()) return false;

  fastjet::PseudoJet lead = jets[0];
  fastjet::JetDefinition subDef(fastjet::antikt_algorithm, R_sub);
  fastjet::ClusterSequence csSub(lead.constituents(), subDef);
  auto subjets = fastjet::sorted_by_pt(csSub.inclusive_jets());
  if (subjets.size() < 2) return false;

  const double p1 = subjets[0].perp();
  const double p2 = subjets[1].perp();
  if (p1 < ptCutSubjet || p2 < ptCutSubjet) return false;

  const double sum = p1 + p2;
  if (sum <= 0.0 || !std::isfinite(sum)) return false;

  zsj   = std::min(p1, p2) / sum;       // 0 ≤ z ≤ 0.5
  theta = subjets[0].delta_R(subjets[1]);
  return std::isfinite(zsj) && std::isfinite(theta);
}
