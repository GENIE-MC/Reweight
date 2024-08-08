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
  size_t nch{};
  for (int i{}; i < event.GetEntries(); ++i) {
    auto particle = event.Particle(i);
    if (particle->Status() == kIStHadronInTheNucleus && particle->Charge()) {
      nch++;
    }
  }
  ret.push_back(nch);

  bool cc = event.Summary()->ProcInfo().IsWeakCC();
  ret.push_back(cc);

  int probe = event.Probe()->Pdg();
  ret.push_back(probe);
  
  int nucleus = event.TargetNucleus()->Pdg();
  ret.push_back(nucleus);


  return ret;
}

bool ObservableHadronization::IsHandled(const EventRecord &event) const {
  return event.FinalStateHadronicSystem(); // for events with HadronicSystem
}

void ObservableHadronization::LoadConfig(void) {}

} // namespace rew
} // namespace genie