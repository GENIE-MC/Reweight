#include "ProfSpline/ObservableSplines.h"
#include "ProfSpline/Ipol.h"
#include <cstddef>

namespace genie {
namespace rew {

ObservableI *ObservableSplines::GetObservable() const {
  return observable.get();
}

const Professor::Ipol &ObservableSplines::GetBin(size_t bin_id) const {
  return bins[bin_id];
}

const Professor::Ipol &ObservableSplines::operator[](size_t bin_id) const {
  return bins[bin_id];
}

const Professor::Ipol &
ObservableSplines::operator[](const std::vector<double> &vars) const {
  return bins[GetObservablesBinID(vars)];
}

double ObservableSplines::GetDXsec(const EventRecord *evt,
                                    const std::vector<double> &para) const {
  return (*this)[observable->HandleEventRecord(evt)].value(para);
}

void ObservableSplines::InitializeBins(
    const std::vector<std::vector<double>> &m_bin_edges) {
  bin_edges = m_bin_edges;
}

size_t ObservableSplines::GetObservablesBinID(
    const std::vector<double> &values) const {
  assert(values.size() == bin_edges.size()); // other wise we are having trouble
  size_t bin_id = 0;
  for (size_t i = 0; i < values.size();
       ++i) { // iterate over all the dimensions
              // depends on how the linear index is calculated we may need to
              // reverse the order here
    bin_id *= bin_edges[i].size() - 1;
    auto value = values[i];
    auto &edges = bin_edges[i];
    bool overflow{true};

    for (size_t j = 0; j < edges.size() - 1; ++j) {
      if (value >= edges[j] && value < edges[j + 1]) {
        bin_id += j;
        overflow = false;
        break;
      }
    }
    if (overflow) {
      bin_id += edges.size() - 2;
      // TODO: do we need to raise an error here?
    }
  }
  return bin_id;
}

void ObservableSplines::InitialIpols(const std::vector<std::string> &lines) {
  bins.clear();
  bins.reserve(lines.size());
  for (const auto &line : lines) {
    bins.emplace_back(line);
  }
}

} // namespace rew
} // namespace genie