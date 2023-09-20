#include "RwCalculators/GReWeightProfessor.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/XmlParserUtils.h"
#include <cstdlib>
#include <fstream>
#include <memory>
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

void GReWeightProfessor::Initialize(std::string conf_file) {
  // This is a very simliar function to the one in
  // GENIE_COMPARISONS/src/Observables/GeneralReweightObs.cxx
  // we may consider merging them in the future
  ReadComparionConf(conf_file);
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
    } else {
      LOG("GReWeightProfessor::ReadProf2Spline", pINFO)
          << "Unknown line: " << line;
    }
  }
  observable_splines->InitializeIpols(var_lines);
}

void GReWeightProfessor::ReadComparionConf(std::string conf_file) {
  xmlDocPtr xmldoc = xmlParseFile(conf_file.c_str());
  if (xmldoc == NULL) {
    LOG("GReWeightProfessor::Initialize", pERROR)
        << "Error parsing file " << conf_file << "\n";
    exit(1);
  }
  dimension = utils::xml::GetInt(xmldoc, "config/observables/dimension");
  std::vector<std::vector<double>> bin_edges{};
  bin_edges.resize(dimension);
  // we want to reference the bin edge file relative to the GENIE_COMPARISONS
  const char *basedir = std::getenv("GENIE_COMPARISONS");
  basedir = basedir ? basedir : "";
  size_t total_bins = 1;
  for (size_t i = 0; i < dimension; ++i) {
    auto nodepath =
        "config/observables/bin_edges/bin_edge_" + std::to_string(i + 1);
    std::string bin_edge_filepath = utils::xml::GetString(xmldoc, nodepath);
    bin_edge_filepath = utils::str::TrimSpaces(bin_edge_filepath);
    auto filepath = (bin_edge_filepath[0] != '/')
                        ? (std::string{basedir} + "/" + bin_edge_filepath)
                        : bin_edge_filepath;
    // LOG("", pINFO) << "reading from " << filepath;
    
    LOG("GReWeightProfessor", pINFO)
        << "Reading bin edges from " << filepath;
    std::ifstream bin_edge_file(filepath);

    std::string line;
    while (std::getline(bin_edge_file, line)) {
      bin_edges[i].push_back(std::stod(line));
    }
    total_bins *= bin_edges[i].size() - 1;
  }
  auto obsname = genie::utils::xml::GetString(xmldoc, "config/observable_name");
  auto disc_bin_names =
      genie::utils::xml::GetString(xmldoc, "config/enabled_discrete_bins");
  auto disc_bin_names_vec = utils::str::Split(disc_bin_names, ",");
  observable_splines = std::make_unique<rew::ObservableSplines>();
  observable_splines->InitializeBins(bin_edges);
  observable_splines->InitializeObservable(obsname);
  observable_splines->InitializeDiscreteBins(disc_bin_names_vec);
  total_bins *= observable_splines->GetNChannel();
  LOG("GReWeightProfessor", pINFO) << "Initized with " << total_bins << "bins";
}

GReWeightProfessor::GReWeightProfessor(std::string name)
    : GReWeightModel(name) {}

} // namespace rew
} // namespace genie
