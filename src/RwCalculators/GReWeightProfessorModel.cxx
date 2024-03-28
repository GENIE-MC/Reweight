#include "RwCalculators/GReWeightProfessorModel.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/XmlParserUtils.h"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

namespace genie {
namespace rew {

//! does the current weight calculator handle the input nuisance param?
bool GReWeightProfessorModel::IsHandled(GSyst_t syst) const { return true; }

//! update the value for the specified nuisance param
void GReWeightProfessorModel::SetSystematic(GSyst_t syst, double val) {
  // TODO: figure out how to do it
  // We need a map of GSyst_t to the index of the nuisance parameter
  // Again this should match whatever in the comparison package
  // And same information will be used in GReWeightProfessor::IsHandled
  // systematics_values[systematics_map[syst]] = val;
}

bool GReWeightProfessorModel::AppliesTo(const EventRecord &event) const {
  return observable->IsHandled(event);
}

//!  set all nuisance parameters to default values
void GReWeightProfessorModel::Reset(void) {}

//! propagate updated nuisance parameter values to actual MC, etc
void GReWeightProfessorModel::Reconfigure(void) {}

//! calculate a weight for the input event using the current nuisance param
//! values
double GReWeightProfessorModel::CalcWeight(const genie::EventRecord &event) {
  // return observable_splines->GetRatio(event, systematics_values, orig_value);
  if (!observable) {
    LOG("GReWeightProfessor::CalcWeight", pERROR)
        << "Cannot find observable splines for event. "
           "Missing Check GReWeightProfessor::AppliesTo ? "
           "But Still return weight 1";

    return 1;
  }
  return observable->GetRatio(event, systematics_values, orig_value);
}

// void GReWeightProfessor::Initialize(std::string conf_file) {
//   // This is a very simliar function to the one in
//   // GENIE_COMPARISONS/src/Observables/GeneralReweightObs.cxx
//   // we may consider merging them in the future
//   // ReadComparionConf(conf_file);
//   // ReadComparionXML(conf_file);
// }

GReWeightProfessorModel::GReWeightProfessorModel(std::string name)
    : GReWeightModel(name) {}

void GReWeightProfessorModel::InitializeIpols(std::vector<std::string> lines) {
  observable->InitializeIpols(lines);
}

} // namespace rew
} // namespace genie
