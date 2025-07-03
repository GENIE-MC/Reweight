#include "ProfSpline/ObservableSplines.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgId.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/XmlParserUtils.h"
#include "Professor/Ipol.h"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>

namespace genie {
namespace rew {

const Professor::Ipol &ObservableSplines::GetBin(size_t bin_id) const {
  return bins[bin_id];
}

double ObservableSplines::GetDXsec(const EventRecord &evt,
                                   const std::vector<double> &para) const {
  auto bin_id = lookupBinID(observable->CalcKinematicVariables(evt));
  return bins[bin_id].value(para);
}

double ObservableSplines::GetRatio(const EventRecord &evt,
                                   const std::vector<double> &para,
                                   const std::vector<double> &para_orig) const {
  if (!observable->IsHandled(evt)) {
    return 1.;
  }
  if (!(observable->ChannelID(evt) == channel)) {
    return 1.;
  }
  auto obvs = observable->CalcKinematicVariables(evt);
  auto new_weight = GetValueInterpolated(obvs, para);
  auto old_weight = GetValueInterpolated(obvs, para_orig);
  const auto product = new_weight / old_weight;
  /// \uml{skip}
  if (product < 0 || isnan(product)) {
    auto bin_id = lookupBinID(obvs);
    LOG("ObservableSplines", pERROR) << "Negative ratio: " << new_weight << " \
    / " << old_weight << " for bin " << bin_id
                                     << " with kinematic variables: " <<
        [](std::vector<double> &vec) {
          std::stringstream ss;
          for (auto i : vec) {
            ss << i << "\t";
          }
          return ss.str();
        }(obvs) << " new weight: " << new_weight
                                     << " old weight: " << old_weight << " \
    ratio: " << new_weight / old_weight;
    return 1;
  }
  return product;
}

void ObservableSplines::InitializeIpols(const std::vector<std::string> &lines) {
  bins.clear();
  bins.reserve(lines.size());
  LOG("ObservableSplines", pNOTICE)
      << "Initializing " << lines.size() << " bins";
  for (const auto &line : lines) {
    if (line.empty()) {
      LOG("ObservableSplines", pERROR)
          << "Empty line in spline definition, this could lead to errors";
      exit(1);
    }
    bins.emplace_back(line);
  }
}

void ObservableSplines::InitializeObservable(const std::string name,
                                             const std::string config) {
  observable = dynamic_cast<const genie::rew::RwgKineSpace *>(
      AlgFactory::Instance()->GetAlgorithm(name, config));
  if (!observable) {
    LOG("ObservableSplines", pFATAL)
        << "Cannot find observable " << name << " in rew algorithm list";
    exit(1);
  }
}

void ObservableSplines::InitializeObservable(const std::string AlgID) {
  auto div = AlgID.find_first_of("/");
  if (div == std::string::npos) {
    LOG("ObservableSplines", pINFO)
        << "Cannot find config for observable " << AlgID << ", using NoConfig";
    InitializeObservable(AlgID, "NoConfig");
    return;
  }

  auto name = AlgID.substr(0, div);
  auto config = AlgID.substr(div + 1);
  LOG("ObservableSplines", pINFO)
      << "Initializing observable " << name << " with config " << config;
  InitializeObservable(name, config);
}

void ObservableSplines::InitializeObservable(AlgId id) {
  observable = dynamic_cast<const genie::rew::RwgKineSpace *>(
      AlgFactory::Instance()->GetAlgorithm(id));
  if (!observable) {
    LOG("ObservableSplines", pFATAL)
        << "Cannot find observable " << id << " in rew algorithm list";
    exit(1);
  }
}

double ObservableSplines::GetValueInterpolated(
    const std::vector<double> &obvs, const std::vector<double> &paras) const {
  auto bin_id = lookupBinID(obvs);
  return bins[bin_id].value(paras); // noop for now
}

size_t ObservableSplines::lookupBinID(const std::vector<double> &obvs) const {
  assert(obvs.size() == dimension);
  for (size_t binid{}; binid < bin_edges.size(); ++binid) {
    const auto bin = bin_edges[binid];
    bool match{true};
    for (size_t i = 0; i < dimension; ++i) {
      if (obvs[i] < bin[i].first || obvs[i] > bin[i].second) {
        match = false;
      }
    }
    if (match)
      return binid;
  }
  throw std::runtime_error("Cannot find bin for observable");
}

} // namespace rew
} // namespace genie