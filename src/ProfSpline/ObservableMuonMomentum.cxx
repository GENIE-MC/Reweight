#include "ProfSpline/ObservableMuonMomentum.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"

namespace genie {
namespace rew {

ObservableMuonMomentum::ObservableMuonMomentum()
    : ObservableI("genie::rew::ObservableMuonMomentum") {}

std::vector<double>
ObservableMuonMomentum::GetKinematicVariables(const EventRecord &event) const {
  std::vector<double> muon_momentum;
  double muon_momentum_v = event.FinalStatePrimaryLepton()->P4()->P();
  muon_momentum.resize(1);
  muon_momentum[0] = muon_momentum_v;
  return muon_momentum;
}

bool ObservableMuonMomentum::IsHandled(const EventRecord &event) const {
  if (event.Summary()->ProcInfo().IsWeakCC() != isCC) {
    return false;
  }
  if (channelid == "ALL") {
    return true;
  } else if (channelid == "QEL") {
    return event.Summary()->ProcInfo().IsQuasiElastic();
  } else if (channelid == "DIS") {
    return event.Summary()->ProcInfo().IsDeepInelastic();
  } else if (channelid == "MEC") {
    return event.Summary()->ProcInfo().IsMEC();
  } else if (channelid == "RES") {
    return event.Summary()->ProcInfo().IsResonant();
  }
  return false;
}

void ObservableMuonMomentum::Configure(const Registry &config) {
  Algorithm::Configure(config);
  LoadConfig();
}
//____________________________________________________________________________
void ObservableMuonMomentum::Configure(string param_set) {
  Algorithm::Configure(param_set);
  LoadConfig();
}

void ObservableMuonMomentum::LoadConfig(void) {
  GetParamDef("channelid", channelid, std::string("ALL"));
  GetParamDef("isCC", isCC, true);
  LOG("ObservableMuonMomentum", pINFO)
      << "Configured with channelid: " << channelid << " isCC: " << isCC;
}

bool ObservableMuonMomentum::IsCC() const { return isCC; }

} // namespace rew
} // namespace genie