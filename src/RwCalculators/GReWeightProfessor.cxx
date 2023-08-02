#include "RwCalculators/GReWeightProfessor.h"
#include <fstream>
#include <sstream>
#include <string>

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
  return observable_splines->GetRatio(event, systematics_values, orig_value);
}

void GReWeightProfessor::Initialize(std::string spline_filepath) {
  observable_splines = std::make_unique<ObservableSplines>();
  ReadProf2Spline(spline_filepath);
  std::string observable_name = "muon_momentum";
  observable_splines->InitializeObservable(observable_name);
  std::vector<std::vector<double>> binning{};
  observable_splines->InitializeBins(binning);
}

void GReWeightProfessor::ReadProf2Spline(std::string filepath) {
  std::ifstream spline_file{filepath};
  std::vector<std::string> var_lines{};
  for (std::string line; std::getline(spline_file, line);) {
    auto seperator = line.find(":");
    if (seperator != std::string::npos) {
      auto name = line.substr(0, seperator);
      auto var = line.substr(seperator + 1);
      if (name == "ParamNames") {
        std::stringstream ss{var};
        for (std::string param; ss >> param;) {
          spline_vars.push_back(param);
        }
      } else if (name == "MinParamVals") {
        std::stringstream ss{var};
        for (double param; ss >> param;) {
          var_min.push_back(param);
        }
      } else if (name == "MaxParamVals") {
        std::stringstream ss{var};
        for (double param; ss >> param;) {
          var_max.push_back(param);
        }
      } else if (name == "Dimension") {
        std::stringstream ss{var};
        ss >> dimension;
      }
    } else if (line.find("#") != std::string::npos) {
      std::string errline{}; // not used now
      auto &varline = var_lines.emplace_back();
      std::getline(spline_file, varline);
      std::getline(spline_file, errline);
    }
  }
  observable_splines->InitializeIpols(var_lines);
}

void GReWeightProfessor::ReadComparionConf(std::string filepath) {}

} // namespace rew
} // namespace genie
