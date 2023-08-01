#include "ProfSpline/ObservableMuonMomentum.h"
#include "Framework/GHEP/GHepParticle.h"
#include "ProfSpline/ObservableF.h"

namespace genie {
namespace rew {

REGISTER_OBSERVABLE("muon_momentum", ObservableMuonMomentum)

std::vector<double>
ObservableMuonMomentum::HandleEventRecord(const EventRecord *event) {
  std::vector<double> muon_momentum;
  double muon_momentum_v = event->FinalStatePrimaryLepton()->P4()->P();
  muon_momentum.resize(1);
  muon_momentum[0] = muon_momentum_v;
  return muon_momentum;
}

} // namespace rew
} // namespace genie