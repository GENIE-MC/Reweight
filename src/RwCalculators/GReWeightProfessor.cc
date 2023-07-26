#include "RwCalculators/GReWeightProfessor.h"

namespace genie {
namespace rew {
bool GReWeightProfessor::AppliesTo(ScatteringType_t type, bool is_cc) const {
  return true;
}

//! does the current weight calculator handle the input nuisance param?
bool GReWeightProfessor::IsHandled(GSyst_t syst) const { return true; }

//! update the value for the specified nuisance param
void GReWeightProfessor::SetSystematic(GSyst_t syst, double val) {
  // TODO: figure out how to do it
  // We need a map of GSyst_t to the index of the nuisance parameter
  // Again this should match whatever in the comparison package
  // And same information will be used in GReWeightProfessor::IsHandled
  // systematics_values[systematics_map[syst]] = val;
}

//!  set all nuisance parameters to default values
void GReWeightProfessor::Reset(void) {}

//! propagate updated nuisance parameter values to actual MC, etc
void GReWeightProfessor::Reconfigure(void) {}

//! calculate a weight for the input event using the current nuisance param
//! values
double GReWeightProfessor::CalcWeight(const genie::EventRecord &event) {
  return observable_splines->get_value(&event, systematics_values);
}
} // namespace rew
} // namespace genie
