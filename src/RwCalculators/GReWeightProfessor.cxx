#include "RwCalculators/GReWeightProfessor.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/XmlParserUtils.h"
#include "ProfSpline/KinematicVariables.h"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>

namespace genie {
namespace rew {

std::tuple<int, int> GetProbTarget(const EventRecord &event) {
  GHepParticle *probe = event.Probe();
  assert(probe);
  Interaction *interaction = event.Summary();
  Target *target = (interaction->InitState()).TgtPtr();
  assert(target);
  auto probid = probe->Pdg();
  auto targetid = target->Pdg();
  return {probid, targetid};
}

bool GReWeightProfessor::AppliesTo(const EventRecord &event) const {
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
  // return observable_splines->GetRatio(event, systematics_values, orig_value);
  double weight{1.};
  for (auto &&calc : observables) {
    weight *= calc.GetRatio(event, orig_value, systematics_values);
  }
  return weight;
}

std::map<std::tuple<std::string /*observable id*/,
                    ChannelIDs /*Channel selection ID*/>,
         std::vector<std::string> /*the vars lines from prof2*/>
GReWeightProfessor::ReadProf2Spline(std::string filepath) {
  std::map<std::tuple<std::string, ChannelIDs>, std::vector<std::string>> ret;
  std::ifstream spline_file{filepath};
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
      // "/${orbname}/${orbname}__<channel_string>#..."
      auto path = line.substr(0, line.find("#"));
      auto name = path.substr(path.find_last_of("/") + 1);
      auto orbname = path.substr(1, path.find_last_of("/") - 1);
      auto sep_loc = name.find_last_of("__");
      auto channelstr = name.substr(sep_loc + 2);
      ChannelIDs channel(channelstr, "_");
      // auto name_trim_head = name.substr(orbname.length() + 1);
      // auto seperator_loc = name_trim_head.find_last_of("_");
      // auto nuclid = std::stoi(name_trim_head.substr(seperator_loc + 1));
      // auto flavor = std::stoi((name_trim_head.substr(0, seperator_loc)));

      std::string errline{}; // not used now
      // auto &varline = var_lines.emplace_back();
      auto &varline = ret[std::make_tuple(orbname, channel)].emplace_back();
      std::getline(spline_file, varline);
      std::getline(spline_file, errline);
    } else {
      LOG("GReWeightProfessor::ReadProf2Spline", pINFO)
          << "Unknown line: " << line;
    }
  }
  return ret;
}

GReWeightProfessor::GReWeightProfessor(std::string name)
    : GReWeightModel(name) {}

void GReWeightProfessor::ReadComparisonXML(std::string filepath,
                                           std::string spline_path) {
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
  auto splines = ReadProf2Spline(spline_path);
  // get node "binning"
  auto node = utils::xml::FindNode(doc, "binning");
  for (auto observable_node = node->children; observable_node;
       observable_node = observable_node->next) {
    if (observable_node->type == XML_ELEMENT_NODE) {
      auto nodename = std::string((const char *)observable_node->name);
      auto algid = utils::xml::GetAttribute(observable_node, "Algid");
      for (auto blocknode = observable_node->children; blocknode;
           blocknode = blocknode->next) {
        if (blocknode->type != XML_ELEMENT_NODE) {
          continue;
        }
        ChannelIDs channel{};
        // TODO: read channel from file
        auto name = utils::xml::GetAttribute(blocknode, "name");
        // auto prob = std::stoi(utils::xml::GetAttribute(blocknode, "prob"));
        // auto nuclid = std::stoi(utils::xml::GetAttribute(blocknode, "nucl"));
        auto bin_count =
            std::stoul(utils::xml::GetAttribute(blocknode, "size"));
        // auto bin_count =
        //     std::stoul(utils::xml::GetAttribute(blocknode, "size"));
        auto dimension =
            std::stoul(utils::xml::GetAttribute(blocknode, "dimension"));
        std::vector<std::vector<std::pair<double, double>>> bin_edges{};
        bin_edges.resize(bin_count);
        std::vector<std::set<size_t>> first_neighbour{};
        first_neighbour.resize(bin_count);
        for (auto cur = blocknode->children; cur; cur = cur->next) {
          if (cur->type != XML_ELEMENT_NODE)
            continue;
          auto bin_id = std::stoul(utils::xml::GetAttribute(cur, "binid"));

          // for each bin iteriate over all the attributes
          // and get the bin edges and first neighbours
          for (auto element = cur->children; element; element = element->next) {
            if (element->type != XML_ELEMENT_NODE)
              continue;
            if (!xmlStrcmp(element->name, (const xmlChar *)"axis")) {
              bin_edges[bin_id].resize(dimension);

              auto axis_id =
                  std::stoul(utils::xml::GetAttribute(element, "axisid"));
              for (auto axis = element->children; axis; axis = axis->next) {
                // std::pair<double, double> axis_range;
                auto &axis_range = bin_edges[bin_id][axis_id];
                if (!xmlStrcmp(axis->name, (const xmlChar *)"min")) {
                  auto str =
                      xmlNodeListGetString(doc, axis->xmlChildrenNode, 1);
                  axis_range.first = std::stod((const char *)str);
                } else if (!xmlStrcmp(axis->name, (const xmlChar *)"max")) {
                  auto str =
                      xmlNodeListGetString(doc, axis->xmlChildrenNode, 1);
                  axis_range.second = std::stod((const char *)str);
                }
                // bin_edges[bin_id][axis_id] = axis_range;
              }
            } else if (!xmlStrcmp(element->name, (const xmlChar *)"neighbor")) {
              auto str = xmlNodeListGetString(doc, element->xmlChildrenNode, 1);
              std::stringstream str_view((const char *)str);
              std::string item;
              while (std::getline(str_view, item, ',')) {
                first_neighbour[bin_id].insert(std::stoul(item));
              }
            }
          }
        } // end of per block loop
        auto &&observable = observables.emplace_back(bin_edges, first_neighbour,
                                                     channel, dimension);
        observable.InitializeObservable(algid);
        observable.InitializeIpols(splines.at(std::make_tuple(name, channel)));
      }
    }
  }
}

} // namespace rew
} // namespace genie
