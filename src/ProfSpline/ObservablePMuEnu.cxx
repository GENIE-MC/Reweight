#include "ProfSpline/ObservablePMuEnu.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"

namespace genie {
namespace rew {
std::vector<double> ObservablePMuEnu::GetKinematicVariables(const EventRecord &event) const {
  // std::vector<double> ret;
  double pmu = event.FinalStatePrimaryLepton()->P4()->P();
  double enu = event.Probe()->E();
  return {pmu, enu};
}

} // namespace rew
} // namespace genie