#include <algorithm>

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "ObservableHadronization.h"
#include "TLorentzVector.h"

namespace genie {
namespace rew {

ObservableHadronization::ObservableHadronization()
    : RwgKineSpace("genie::rew::ObservableHadronization") {}

ObservableHadronization::ObservableHadronization(std::string config)
    : RwgKineSpace("genie::rew::ObservableHadronization", config) {}

KinematicVariables ObservableHadronization::CalcKinematicVariables(
    const EventRecord &event) const {
  std::vector<double> ret;

  double W = event.Summary()->Kine().W(true);
  ret.push_back(W);

  auto y = event.Summary()->Kine().y(true);
  ret.push_back(y);

  return ret;
}

ChannelIDs ObservableHadronization::ChannelID(const EventRecord &event) const {
  ChannelIDs ret{};

  if (event.TargetNucleus()) {
    ret.push_back(event.Probe()->Pdg());
    ret.push_back(event.TargetNucleus()->Pdg());
  } else if (event.HitNucleon()) {
    ret.push_back(event.Probe()->Pdg());
    ret.push_back(event.HitNucleon()->Pdg());
  } else {
    ret.push_back(event.Probe()->Pdg());
    ret.push_back(0);
  }

  const auto nucleus = ret[1];

  size_t nch{}, nn{};
  auto ne = event.GetEntriesFast();
  for (int i{}; i < ne; ++i) {
    auto particle = event.Particle(i);
    if (particle->Status() == ((nucleus == 2212 || nucleus == 2112)
                                   ? kIStStableFinalState
                                   : kIStHadronInTheNucleus)) {
      nch += (particle->Charge() != 0);
      nn += (particle->Charge() == 0);
    }
  }
  nch = std::min<size_t>(nch, 17);
  ret.push_back(nch);

  nn = std::min<size_t>(nn, 3);
  ret.push_back(nn);

  bool cc = event.Summary()->ProcInfo().IsWeakCC();
  ret.push_back(cc);

  // int probe = event.Probe()->Pdg();
  // ret.push_back(probe);

  // int nucleus = event.TargetNucleus()->Pdg();
  // ret.push_back(nucleus);

  return ret;
}

bool ObservableHadronization::IsHandled(const EventRecord &event) const {
  return event.Summary()->ProcInfo().IsResonant() ||
         event.Summary()->ProcInfo().IsDeepInelastic();
}

void ObservableHadronization::LoadConfig(void) {}

} // namespace rew
} // namespace genie