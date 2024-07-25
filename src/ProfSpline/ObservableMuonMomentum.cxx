#include "ProfSpline/ObservableMuonMomentum.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"

namespace genie {
namespace rew {

ObservableMuonMomentum::ObservableMuonMomentum()
    : RwgKineSpace("genie::rew::ObservableMuonMomentum") {}

ObservableMuonMomentum::ObservableMuonMomentum(std::string config)
    : RwgKineSpace("genie::rew::ObservableMuonMomentum", config) {}

KinematicVariables
ObservableMuonMomentum::CalcKinematicVariables(const EventRecord &event) const {
  // std::vector<double> muon_momentum;
  double muon_momentum_v = event.FinalStatePrimaryLepton()->P4()->P();
  return {muon_momentum_v};
}

ChannelIDs ObservableMuonMomentum::ChannelID(const EventRecord &event) const {
  if (event.TargetNucleus()) {
    return {
        event.Probe()->Pdg(),
        event.TargetNucleus()->Pdg(),
    };
  } else if (event.HitNucleon()) {
    return {
        event.Probe()->Pdg(),
        event.HitNucleon()->Pdg(),
    };
  } else {
    return {
        event.Probe()->Pdg(),
        0,
    };
  }
}

bool ObservableMuonMomentum::IsHandled(const EventRecord &event) const {
  if (!event.Summary()->ProcInfo().IsWeakCC()) {
    return false;
  }
  if (channelid == "ALL") {
    return true;
  }
  return ScatteringType::AsString(
             event.Summary()->ProcInfo().ScatteringTypeId()) == channelid;
}

void ObservableMuonMomentum::LoadConfig(void) {
  GetParamDef("channelid", channelid, std::string("ALL"));
}

} // namespace rew
} // namespace genie