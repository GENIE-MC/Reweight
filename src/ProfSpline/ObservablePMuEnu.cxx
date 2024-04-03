#include "ProfSpline/ObservablePMuEnu.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"

namespace genie {
namespace rew {
ObservablePMuEnu::ObservablePMuEnu()
    : RwgKineSpace("genie::rew::ObservablePMuEnu") {}

ObservablePMuEnu::ObservablePMuEnu(std::string config)
    : RwgKineSpace("genie::rew::ObservablePMuEnu", config) {}

KinematicVariables
ObservablePMuEnu::CalcKinematicVariables(const EventRecord &event) const {
  KinematicVariables ret{};
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
  auto &vars = ret.GetVars();
  vars.resize(2);
  vars[0] = pmu;
  vars[1] = enu;

  auto &channel = ret.GetChannel();
  channel.resize(2);
  int nucleon = event.TargetNucleus()->Pdg();
  int neutrino = event.Probe()->Pdg();
  channel[0] = nucleon;
  channel[1] = neutrino;
  return ret;
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