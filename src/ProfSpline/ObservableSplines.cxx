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
  // // LOG("ObservableSplines", pINFO) << "Bin ID: " << bin_id;
  // if (bin_id >= bins.size()) {
  //   LOG("ObservableSplines", pERROR) << "Bin ID out of range: " << bin_id;
  //   return 1;
  // }
  // auto &&bin = GetBin(bin_id);
  // auto new_weight = bin.value(para);
  // auto old_weight = bin.value(para_orig);
  auto channelid = GetChannelID(evt);
  // auto bin_id = GetObservablesBinID(evt);
  auto obvs = observable->GetKinematicVariables(evt);
  auto new_weight = GetValueInterpolated(channelid, obvs, para);
  auto old_weight = GetValueInterpolated(channelid, obvs, para_orig);
  const auto product = new_weight / old_weight;
  if (product < 0 || isnan(product)) {
    // Sometimes we can reach here, this can be due to rare
    // events that not being recorded during generation of
    // splines
    // OR
    // Something went wrong during generation of spline or doing
    // intepolation
    // .. For now we just "passthough" the events
    // .. But warn here since too much of such events
    // .. is definitely a problem
    auto &&kin_vars = observable->GetKinematicVariables(evt);
    LOG("ObservableSplines", pERROR) << "Negative ratio: " << new_weight << " \
    / " << old_weight << " for bin " << bin_id
                                     << " with kinematic variables: " <<
        [](std::vector<double> &vec) {
          std::stringstream ss;
          for (auto i : vec) {
            ss << i << "\t";
          }
          return ss.str();
        }(kin_vars)
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

double ObservableSplines::GetValueInterpolated(
    size_t channel_id, const std::vector<double> &obvs,
    const std::vector<double> &paras) const {
  auto var_loc = binning.GetObservablesBinLoc(obvs);
  auto center_bin_id = binning.GetObservablesBinID(obvs);
  const auto center_value =
      bins[toBinID(channel_id, center_bin_id)].value(paras);
  double diffsum{};
  for (size_t i = 0; i < var_loc.size(); ++i) {
    auto near_bin_array = center_bin_id;
    auto nbins = binning.GetAxis(i).GetNbins();
    auto bin_id = binning.GetAxis(i).FindBin(obvs[i]);
    auto bin_center = binning.GetAxis(i).GetBinCenter(bin_id);
    int near_bin_id{};
    if (bin_id == 1) {
      near_bin_id = 2;
    } else if (bin_id == nbins) {
      near_bin_id = nbins - 1;
    } else {
      if ((abs(obvs[i] - binning.GetAxis(i).GetBinCenter(bin_id - 1)) >
           abs(obvs[i] - binning.GetAxis(i).GetBinCenter(bin_id + 1)))) {
        near_bin_id = bin_id + 1;
      } else {
        near_bin_id = bin_id - 1;
      }
    }
    auto near_bin_center = binning.GetAxis(i).GetBinCenter(near_bin_id);
    near_bin_array[i] = near_bin_id;
    auto near_bin_value =
        bins[toBinID(channel_id, near_bin_array)].value(paras);
    auto diff = (near_bin_value - center_value) /
                (near_bin_center - bin_center) * (obvs[i] - bin_center);
    diffsum += diff;
  }
  return center_value + diffsum;
}

size_t ObservableSplines::toBinID(size_t channel_id,
                                  std::vector<int> bin_ids) const {
  return binning.GetObservablesBinIDLinearized(bin_ids) * GetNChannel() +
         channel_id;
}

} // namespace rew
} // namespace genie