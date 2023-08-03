#include "ProfSpline/ObservableMuonMomentum.h"
#include "Framework/GHEP/GHepParticle.h"

namespace genie {
namespace rew {


std::vector<double>
ObservableMuonMomentum::GetKinematicVariables(const EventRecord &event) const {
  std::vector<double> muon_momentum;
  double muon_momentum_v = event.FinalStatePrimaryLepton()->P4()->P();
  muon_momentum.resize(1);
  muon_momentum[0] = muon_momentum_v;
  return muon_momentum;
}

} // namespace rew
} // namespace genie