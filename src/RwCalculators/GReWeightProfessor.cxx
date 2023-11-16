#include "RwCalculators/GReWeightProfessor.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/XmlParserUtils.h"
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

namespace genie {
namespace rew {
bool GReWeightProfessor::AppliesTo(const EventRecord & event) const {
  bool ret {false};
  for (auto &&[id, obs] : observable_map_from_id) {
    // if 
  }
  return ret;
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
  // return observable_splines->GetRatio(event, systematics_values, orig_value);
}

void GReWeightProfessor::Initialize(std::string conf_file) {
  // This is a very simliar function to the one in
  // GENIE_COMPARISONS/src/Observables/GeneralReweightObs.cxx
  // we may consider merging them in the future
  // ReadComparionConf(conf_file);
}

void GReWeightProfessor::ReadProf2Spline(std::string filepath) {
  std::ifstream spline_file{filepath};
  // std::vector<std::string> var_lines{};
  std::map<std::tuple<std::string, int, int>, std::vector<std::string>>
      var_lines{};
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
      }
    } else if (line.find("#") != std::string::npos) {
      // "/${orbname}/${flavor}_${nuclid}#..."
      auto path = line.substr(0, line.find("#"));
      auto name = path.substr(path.find_last_of("/") + 1);
      auto orbname = path.substr(1, path.find_last_of("/"));
      auto seperator_loc = name.find_last_of("_");
      auto nuclid = std::stoi(name.substr(seperator_loc + 1));
      auto flavor = std::stoi((name.substr(0, seperator_loc)));

      std::string errline{}; // not used now
      // auto &varline = var_lines.emplace_back();
      auto &varline =
          var_lines[std::make_tuple(orbname, flavor, nuclid)]
              .emplace_back();
      std::getline(spline_file, varline);
      std::getline(spline_file, errline);
    } else {
      LOG("GReWeightProfessor::ReadProf2Spline", pINFO)
          << "Unknown line: " << line;
    }
  }
  for (auto &&[id, lines] : var_lines) {
    auto &observable = observable_map_from_id[id];
    if (!observable) {
      LOG("GReWeightProfessor::ReadProf2Spline", pFATAL)
          << "Cannot find observable " << std::get<0>(id)
          << "for neutrino flavor" << std::get<1>(id) << " and nuclid "
          << std::get<2>(id) << " in rew algorithm list";
    }
    observable->InitializeIpols(lines);
  }
}

GReWeightProfessor::GReWeightProfessor(std::string name)
    : GReWeightModel(name) {}

void GReWeightProfessor::ReadComparionXML(std::string filepath) {
  auto doc = xmlParseFile(filepath.c_str());
  if (!doc) {
    LOG("GReWeightProfessor::ReadComparionXML", pFATAL)
        << "Cannot parse xml file " << filepath;
    exit(1);
  }
  auto root = xmlDocGetRootElement(doc);
  if (!root) {
    LOG("GReWeightProfessor::ReadComparionXML", pFATAL)
        << "Cannot get root element of xml file " << filepath;
    exit(1);
  }
  // get node "binning"
  auto node = utils::xml::FindNode(doc, "binning");
  for (auto observable_node = node->children; observable_node;
       observable_node = observable_node->next) {
    if (observable_node->type == XML_ELEMENT_NODE) {
      auto nodename = std::string((const char *)observable_node->name);
      auto algid = utils::xml::GetAttribute(observable_node, "AlgID");
      for (auto blocknode = observable_node->children; blocknode;
           blocknode = blocknode->next) {
        auto name = utils::xml::GetAttribute(blocknode, "name");
        auto prob = std::stoi(utils::xml::GetAttribute(blocknode, "prob"));
        auto nuclid = std::stoi(utils::xml::GetAttribute(blocknode, "nucl"));
        // auto bin_count =
        //     std::stoul(utils::xml::GetAttribute(blocknode, "size"));
        auto dimension =
            std::stoul(utils::xml::GetAttribute(blocknode, "dimension"));
        std::vector<std::vector<std::pair<double, double>>> bin_edges{};
        std::vector<std::set<size_t>> first_neighbour{};
        for (auto cur = blocknode->children; cur; cur = cur->next) {
          if (xmlStrcmp(cur->name, (const xmlChar *)"bin")) {
            auto bin_id = std::stoul(utils::xml::GetAttribute(cur, "binid"));

            // for each bin iteriate over all the attributes
            // and get the bin edges and first neighbours
            for (auto element = cur->children; element;
                 element = element->next) {
              if (xmlStrcmp(element->name, (const xmlChar *)"axis")) {
                bin_edges[bin_id].resize(dimension);

                auto axis_id =
                    std::stoul(utils::xml::GetAttribute(element, "axisid"));
                for (auto axis = element->children; axis; axis = axis->next) {
                  std::pair<double, double> axis_range;
                  if (xmlStrcmp(axis->name, (const xmlChar *)"min")) {
                    auto str =
                        xmlNodeListGetString(doc, axis->xmlChildrenNode, 1);
                    axis_range.first = std::stod((const char *)str);
                  } else if (xmlStrcmp(axis->name, (const xmlChar *)"max")) {
                    auto str =
                        xmlNodeListGetString(doc, axis->xmlChildrenNode, 1);
                    axis_range.second = std::stod((const char *)str);
                  }
                  bin_edges[bin_id][axis_id] = axis_range;
                }
              } else if (xmlStrcmp(element->name,
                                   (const xmlChar *)"neighbor")) {
                auto str =
                    xmlNodeListGetString(doc, element->xmlChildrenNode, 1);
                std::stringstream str_view((const char *)str);
                std::string item;
                while (std::getline(str_view, item, ',')) {
                  first_neighbour[bin_id].insert(std::stoul(item));
                }
              }
            }
          }
        } // end of per block loop
        auto obj =
            std::make_shared<ObservableSplines>(bin_edges, first_neighbour);
        // observable_map_from_name[name] = obj;
        observable_map_from_id[std::make_tuple(nodename, prob, nuclid)] = obj;
      }
    }
  }
}

} // namespace rew
} // namespace genie
