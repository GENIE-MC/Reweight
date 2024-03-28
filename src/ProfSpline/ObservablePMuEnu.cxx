#include "ProfSpline/ObservablePMuEnu.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"

namespace genie {
namespace rew {
ObservablePMuEnu::ObservablePMuEnu()
    : RwgKineSpace("genie::rew::ObservablePMuEnu") {}

ObservablePMuEnu::ObservablePMuEnu(std::string config)
    : RwgKineSpace("genie::rew::ObservablePMuEnu", config) {}

std::vector<double>
ObservablePMuEnu::KinematicVariables(const EventRecord &event) const {
  // std::vector<double> ret;
  GHepParticle *target_nucleus_p = event.HitNucleon();
  // in case we are dealing with COH, which comes with no hit nucleon
  // we use the target nucleus instead
  // This is equivalent to not doing any boost
  auto target_nucleus =
      target_nucleus_p ? *target_nucleus_p : *(event.TargetNucleus());
  auto boost_vec = -target_nucleus.P4()->BoostVector();
  auto final_state_lepton = *(event.FinalStatePrimaryLepton()->P4());
  final_state_lepton.Boost(boost_vec);
  double pmu = final_state_lepton.P();
  auto probe = *(event.Probe()->P4());
  probe.Boost(boost_vec);
  double enu = probe.E();
  return {pmu, enu};
}

bool ObservablePMuEnu::IsHandled(const EventRecord &event) const {
  if (!event.Summary()->ProcInfo().IsWeakCC()) {
    return false;
  }
  if (channelid == "ALL") {
    return true;
  }
  return ScatteringType::AsString(
             event.Summary()->ProcInfo().ScatteringTypeId()) == channelid;
}

void ObservablePMuEnu::LoadConfig(void) {
  GetParamDef("channelid", channelid, std::string("ALL"));
}

} // namespace rew
} // namespace genie