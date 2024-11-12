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
  // This is a dummy implementation
  std::vector<double> ret;

  double W = event.Summary()->Kine().W(true);
  ret.push_back(W);

  auto np = event.GetEntriesFast();
  TLorentzVector had_system{};
  for (int i{}; i < np; i++) {
    auto particle = event.Particle(i);
    if (particle->Status() == kIStStableFinalState) {
      had_system += *particle->P4();
    }
  }

  auto had_system_boost = had_system.BoostVector();

  double max_p_pion{};
  double sum_of_transverse_momentum{};
  for (int i{}; i < np; i++) {
    auto particle = event.Particle(i);
    if (particle->Status() == kIStStableFinalState &&
        (particle->Pdg() == 211 || particle->Pdg() == -211 ||
         particle->Pdg() == 111)) {
      // auto p = particle->P4()->P();
      // max_p_pion = std::max(max_p_pion, p);
      auto p_in_had_rest_frame = *(particle->P4());
      p_in_had_rest_frame.Boost(-had_system_boost);
      max_p_pion = std::max(max_p_pion, p_in_had_rest_frame.P());
      sum_of_transverse_momentum += p_in_had_rest_frame.Pt();
    }
  }
  if (max_p_pion == 0) {
    max_p_pion = 1e-6;
  }
  // ret.push_back(max_p_pion);
  ret.push_back(sum_of_transverse_momentum);

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

} // namespace rew
} // namespace genie