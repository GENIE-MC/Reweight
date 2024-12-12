#include <algorithm>

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "ObservableHadronization.h"
#include "TLorentzVector.h"

namespace genie::rew {

// ObservableHadronization::ObservableHadronization()
//     : RwgKineSpace("genie::rew::ObservableHadronization") {}

// ObservableHadronization::ObservableHadronization(std::string config)
//     : RwgKineSpace("genie::rew::ObservableHadronization", config) {}

KinematicVariables ObservableHadronization::CalcKinematicVariables(
    const EventRecord &event) const {
  // This is a dummy implementation
  std::vector<double> ret;

  double W = event.Summary()->Kine().W(true);
  ret.push_back(W);
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

  return ret;
}

bool ObservableHadronization::IsHandled(const EventRecord &event) const {
  if (!(event.Summary()->ProcInfo().IsResonant() ||
        event.Summary()->ProcInfo().IsDeepInelastic())) {
    return false;
  }

  // exclude event with Kaons
  for (int i{}; i < event.GetEntriesFast(); ++i) {
    auto particle = event.Particle(i);
    if (particle->Status() == kIStStableFinalState &&
        pdg::IsKaon(particle->Pdg())) {
      return false;
    }
  }

  return true;
}

void ObservableHadronization::LoadConfig(void) {}

ObservableHadronizationLowW::ObservableHadronizationLowW()
    : ObservableHadronization{"genie::rew::ObservableHadronizationLowW"} {}

ObservableHadronizationLowW::ObservableHadronizationLowW(std::string config)
    : ObservableHadronization{"genie::rew::ObservableHadronizationLowW",
                              config} {}

KinematicVariables ObservableHadronizationLowW::CalcKinematicVariables(
    const EventRecord &event) const {
  auto ret = ObservableHadronization::CalcKinematicVariables(event);
  double visE{};
  auto nentries = event.GetEntriesFast();
  for (int i = 0; i < nentries; i++) {
    auto part = event.Particle(i);
    if (part->Status() == genie::kIStStableFinalState && part->Charge() != 0) {
      // remove mass part for p/n
      visE += part->E();
      if (pdg::IsNeutron(part->Pdg()) || pdg::IsProton(part->Pdg())) {
        visE -= part->Mass();
      }
    }
  }
  ret.push_back(visE);
  return ret;
}

bool ObservableHadronizationLowW::IsHandled(const EventRecord &event) const {
  if (!ObservableHadronization::IsHandled(event)) {
    return false;
  }

  return event.Summary()->Kine().W(true) < 3;
}

ObservableHadronizationHighW::ObservableHadronizationHighW()
    : ObservableHadronization{"genie::rew::ObservableHadronizationHighW"} {}

ObservableHadronizationHighW::ObservableHadronizationHighW(std::string config)
    : ObservableHadronization{"genie::rew::ObservableHadronizationHighW",
                              config} {}

KinematicVariables ObservableHadronizationHighW::CalcKinematicVariables(
    const EventRecord &event) const {
  auto ret = ObservableHadronization::CalcKinematicVariables(event);
  auto np = event.GetEntriesFast();
  TLorentzVector had_system{};
  for (int i{}; i < np; i++) {
    auto particle = event.Particle(i);
    if (particle->Status() == kIStStableFinalState &&
        pdg::IsHadron(particle->Pdg())) {
      had_system += *particle->P4();
    }
  }
  auto had_system_boost = had_system.BoostVector();
  // auto had_system_dir = had_system.Vect().Unit();

  // double sum_of_transverse_momentum{};
  double p_leading_pion{};
  for (int i{}; i < np; i++) {
    auto particle = event.Particle(i);
    if (particle->Status() == kIStStableFinalState &&
        (particle->Pdg() == 211 || particle->Pdg() == -211 ||
         particle->Pdg() == 111)) {
      auto p = *(particle->P4());
      p.Boost(-had_system_boost);
      p_leading_pion = std::max(p_leading_pion, p.P());
    }
  }
  ret.push_back(p_leading_pion);
  return ret;
}

bool ObservableHadronizationHighW::IsHandled(const EventRecord &event) const {
  if (!ObservableHadronization::IsHandled(event)) {
    return false;
  }

  return event.Summary()->Kine().W(true) >= 3;
}

} // namespace genie::rew
