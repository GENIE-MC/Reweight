#include "ObservableHadronization.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Utils/KineUtils.h"

namespace genie {
namespace rew {

ObservableHadronization::ObservableHadronization()
    : RwgKineSpace("genie::rew::ObservableHadronization") {}

ObservableHadronization::ObservableHadronization(std::string config)
    : RwgKineSpace("genie::rew::ObservableHadronization", config) {}

KinematicVariables ObservableHadronization::CalcKinematicVariables(
    const EventRecord &event) const {
  // This is a dummy implementation
  std::vector<double> ret;

  auto interaction = event.Summary();
  auto W = genie::utils::kinematics::W(interaction);
  ret.push_back(W);

  auto hadron = event.FinalStateHadronicSystem();
  auto pt = hadron->P4()->Pt();
  ret.push_back(pt);

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

  size_t nch{};
  for (int i{}; i < event.GetEntries(); ++i) {
    auto particle = event.Particle(i);
    if (particle->Status() == ((nucleus == 2212 || nucleus == 2112)
                                   ? kIStStableFinalState
                                   : kIStHadronInTheNucleus) &&
        particle->Charge()) {
      nch++;
    }
  }
  // if (nch > 4)
  //   nch = 4;
  ret.push_back(nch);

  bool cc = event.Summary()->ProcInfo().IsWeakCC();
  ret.push_back(cc);

  // int probe = event.Probe()->Pdg();
  // ret.push_back(probe);

  // int nucleus = event.TargetNucleus()->Pdg();
  // ret.push_back(nucleus);

  return ret;
}

bool ObservableHadronization::IsHandled(const EventRecord &event) const {
  return event.FinalStateHadronicSystem(); // for events with HadronicSystem
}

void ObservableHadronization::LoadConfig(void) {}

} // namespace rew
} // namespace genie