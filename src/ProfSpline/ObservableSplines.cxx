#include "ProfSpline/ObservableSplines.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Messenger/Messenger.h"
#include "Professor/Ipol.h"
#include <cstddef>

namespace genie {
namespace rew {

const Professor::Ipol &ObservableSplines::GetBin(size_t bin_id) const {
  return bins[bin_id];
}

double ObservableSplines::GetDXsec(const EventRecord &evt,
                                   const std::vector<double> &para) const {
  auto bin_id = GetObservablesBinID(observable->GetKinematicVariables(evt));
  return bins[bin_id].value(para);
}

double ObservableSplines::GetRatio(const EventRecord &evt,
                                   const std::vector<double> &para,
                                   const std::vector<double> &para_orig) const {
  auto bin_id = GetObservablesBinID(observable->GetKinematicVariables(evt));
  return bins[bin_id].value(para) / bins[bin_id].value(para_orig);
}

void ObservableSplines::InitializeBins(
    const std::vector<std::vector<double>> &m_bin_edges) {
  binning.InitializeBins(m_bin_edges);
}

size_t ObservableSplines::GetObservablesBinID(
    const std::vector<double> &values) const {
  return binning.GetObservablesBinIDLinearized(values);
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

} // namespace rew
} // namespace genie