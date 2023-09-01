#include "ProfSpline/ObservableSplines.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Messenger/Messenger.h"
#include "Professor/Ipol.h"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <sstream>

namespace genie {
namespace rew {

const Professor::Ipol &ObservableSplines::GetBin(size_t bin_id) const {
  return bins[bin_id];
}

double ObservableSplines::GetDXsec(const EventRecord &evt,
                                   const std::vector<double> &para) const {
  auto bin_id = GetObservablesBinID(evt);
  return bins[bin_id].value(para);
}

double ObservableSplines::GetRatio(const EventRecord &evt,
                                   const std::vector<double> &para,
                                   const std::vector<double> &para_orig) const {
  // auto bin_id = GetObservablesBinID(kin_vars);
  auto bin_id = GetObservablesBinID(evt);
  // LOG("ObservableSplines", pINFO) << "Bin ID: " << bin_id;
  if (bin_id >= bins.size()) {
    LOG("ObservableSplines", pERROR) << "Bin ID out of range: " << bin_id;
    return 1;
  }
  auto &&bin = GetBin(bin_id);
  auto new_weight = bin.value(para);
  auto old_weight = bin.value(para_orig);
  const auto product = new_weight / old_weight;
  auto vec2str = [](std::vector<double> &vec) {
    std::stringstream ss;
    for (auto i : vec) {
      ss << i << "\t";
    }
    return ss.str();
  };
  if (product < 0 || isnan(product)) {
    auto &&kin_vars = observable->GetKinematicVariables(evt);
    LOG("ObservableSplines", pERROR) << "Negative ratio: " << new_weight << " \
    / " << old_weight << " for bin " << bin_id
                                     << " with kinematic variables: "
                                     << vec2str(kin_vars)
                                     << " new weight: " << new_weight
                                     << " old weight: " << old_weight << " \
    ratio: " << new_weight / old_weight;
    // exit(1);
    return 1;
  }
  return product;
}

void ObservableSplines::InitializeBins(
    const std::vector<std::vector<double>> &m_bin_edges) {
  binning.InitializeBins(m_bin_edges);
}

// size_t ObservableSplines::GetObservablesBinID(
//     const std::vector<double> &values) const {
//   return binning.GetObservablesBinIDLinearized(values);
// }

size_t ObservableSplines::GetObservablesBinID(const EventRecord &event) const {
  return binning.GetObservablesBinIDLinearized(
             observable->GetKinematicVariables(event)) *
             GetNChannel() +
         GetChannelID(event);
}

void ObservableSplines::InitializeIpols(const std::vector<std::string> &lines) {
  bins.clear();
  bins.reserve(lines.size());
  LOG("ObservableSplines", pNOTICE)
      << "Initializing " << lines.size() << " bins";
  for (const auto &line : lines) {
    bins.emplace_back(line);
  }
}

void ObservableSplines::InitializeObservable(const std::string name,
                                             const std::string config) {
  observable = dynamic_cast<const genie::rew::ObservableI *>(
      AlgFactory::Instance()->GetAlgorithm(name, config));
  if (!observable) {
    LOG("ObservableSplines", pFATAL)
        << "Cannot find observable " << name << " in rew algorithm list";
    exit(1);
  }
}

void ObservableSplines::InitializeDiscreteBins(
    const std::vector<std::string> &enabled_bin_names) {
  if (!enabled_bin_names.empty())
    discrete_bins = std::make_unique<ObservableDiscreteBins>(enabled_bin_names);
}

size_t ObservableSplines::GetChannelID(const EventRecord &event) const {
  if (discrete_bins) {
    return discrete_bins->GetBinID(event);
  }
  return 0;
}

size_t ObservableSplines::GetNChannel() const {
  if (discrete_bins) {
    return discrete_bins->GetNBinMax();
  }
  return 1;
}

// size_t ObservableSplines::GetChannelID(const EventRecord &) const {
//   std::array<std::function<bool(const EventRecord &)>, 2> channel{
//       [](const EventRecord &evt) {
//         return evt.Summary()->ProcInfo().IsWeakCC();
//       },
//       [](const EventRecord &evt) {
//         return evt.Summary()->ProcInfo().IsWeakNC();
//       }};

// }

} // namespace rew
} // namespace genie