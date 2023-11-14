#include "ProfSpline/ObservablePMuEnuW.h"
#include "Framework/GHEP/GHepParticle.h"
#include "TVector3.h"

namespace genie {
namespace rew {
std::vector<double>
ObservablePMuEnuW::GetKinematicVariables(const EventRecord &event) const {
  // std::vector<double> ret;
  // GHepParticle *target_nucleus_p = event.HitNucleon();
  // in case we are dealing with COH, which comes with no hit nucleon
  // we use the target nucleus instead
  // This is equivalent to not doing any boost
  // auto target_nucleus =
  //     target_nucleus_p ? *target_nucleus_p : *(event.TargetNucleus());
  // auto boost_vec = -target_nucleus.P4()->BoostVector();
  TVector3 boost_vec(0, 0, 0);
  auto final_state_lepton = *(event.FinalStatePrimaryLepton()->P4());
  final_state_lepton.Boost(boost_vec);
  double pmu = final_state_lepton.P();
  auto probe = *(event.Probe()->P4());
  probe.Boost(boost_vec);
  double enu = probe.E();
  double W = event.Summary()->Kine().W(true);
  if (W == -99999) {
    W = event.Summary()->Kine().W(false);
  }
  W = W == -99999 ? 0 : W;
  return {pmu, enu, W};
}

} // namespace rew
} // namespace genie