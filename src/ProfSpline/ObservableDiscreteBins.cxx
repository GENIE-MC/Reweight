#include "ProfSpline/ObservableDiscreteBins.h"

namespace genie {
namespace rew {
size_t ObservableDiscreteBinCCNC::GetBinID(const EventRecord &evt) const {
  if (evt.Summary()->ProcInfo().IsWeakCC()) {
    return 0;
  }
  if (evt.Summary()->ProcInfo().IsWeakNC()) {
    return 1;
  }
  return 2;
}
size_t ObservableDiscreteBinCCNC::GetNBinMax() const { return 3; }

size_t
ObservableDiscreteBinInteraction::GetBinID(const EventRecord &evt) const {
  if (evt.Summary()->ProcInfo().IsQuasiElastic()) {
    return 0;
  }
  if (evt.Summary()->ProcInfo().IsResonant()) {
    return 1;
  }
  if (evt.Summary()->ProcInfo().IsMEC()) {
    return 2;
  }
  if (evt.Summary()->ProcInfo().IsDeepInelastic()) {
    return 3;
  }
  return 4;
}
size_t ObservableDiscreteBinInteraction::GetNBinMax() const { return 5; }

ObservableDiscreteBins::ObservableDiscreteBins(
    std::vector<std::string> enabled_bin_names) {
  nbin_max = 1;
  bins.reserve(enabled_bin_names.size());
  for (const auto &name : enabled_bin_names) {
    auto alg = dynamic_cast<const ObservableDiscreteBinI *>(
        AlgFactory::Instance()->GetAlgorithm(name));
    if (!alg) {
      LOG("ObservableDiscreteBins", pFATAL)
          << "Cannot find observable " << name << " in rew algorithm list"
          << " as ObservableDiscreteBinI";
      exit(1);
    }
    bins.push_back(alg);
    nbin_max *= alg->GetNBinMax();
  }
}

size_t ObservableDiscreteBins::GetBinID(const EventRecord &evt) const {
  size_t bin_id = 0;
  for (const auto &bin : bins) {
    bin_id *= bin->GetNBinMax();
    bin_id += bin->GetBinID(evt);
  }
  return bin_id;
}

size_t ObservableDiscreteBins::GetNBinMax() const { return nbin_max; }

} // namespace rew
} // namespace genie