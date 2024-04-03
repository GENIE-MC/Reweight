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
  KinematicVariables ret{};
  auto & muon_momentum = ret.GetVars();
  auto & channel = ret.GetChannel();
  double muon_momentum_v = event.FinalStatePrimaryLepton()->P4()->P();
  muon_momentum.resize(1);
  muon_momentum[0] = muon_momentum_v;
  
  channel.resize(2);
  int nucleon = event.TargetNucleus()->Pdg();
  int neutrino = event.Probe()->Pdg();
  channel[0] = nucleon;
  channel[1] = neutrino;
  return ret;
  // return muon_momentum;
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