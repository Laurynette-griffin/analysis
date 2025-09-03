
#ifndef ZSJTOOLS_H
#define ZSJTOOLS_H

#include <vector>
#include <utility>

// Forward decls
class Jet;
class TowerInfoContainer;
class TowerInfo;
class RawTowerGeomContainer;
class TowerBackground;
class PHG4TruthInfoContainer;

namespace fastjet { class PseudoJet; }

// BuildPseudoJets builds FastJet constituents from a Jet's tower components.
// UE modulation: UE_i(eta) * (1 + 2 v2 cos(2(phi-psi2)))
// If doUnsub == true, subtract the modulated UE from tower energy; else leave as-is or add as needed.
std::vector<fastjet::PseudoJet>
BuildPseudoJets(Jet* jet,
                TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
                TowerBackground* bg, float v2, float psi2, bool doUnsub);

// Compute z_sj and theta_sj for a given Jet. Returns true if computed (>=2 subjets with pt>=ptCutSubjet).
// R_jet is the reclustering radius for the jet (e.g., 0.4).
// R_sub is the subjet reclustering radius (e.g., 0.1).
bool ComputeZsjForJet(Jet* jet,
                      TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                      RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
                      TowerBackground* bg, float v2, float psi2,
                      bool doUnsub,
                      double R_jet, double R_sub, double ptCutSubjet,
                      double& zsj, double& theta);

bool ComputeZsjForTruthJet(Jet* truthJet,
                           PHG4TruthInfoContainer* truthInfo,
                           double R_jet, double R_sub, double ptCutSubjet,
                           double& zsj, double& theta);
#endif
