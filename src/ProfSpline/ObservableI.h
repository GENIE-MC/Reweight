// Class to handle all observables as a whole
// This is a interface class

#ifndef _OBSERVABLE_I_
#define _OBSERVABLE_I_
#include "Framework/EventGen/EventRecord.h"
#include <string>
#include <vector>
namespace genie {
namespace rew {
class ObservableI {
public:
  // virtual void Configure(/*some way of configurating*/) = 0;

  // calculate the value for each observable
  // maybe we can merge with binning lookup function below
  // but seperating them may benefit the idea of merging implementation
  // of corresponding ObservablePrediction in comparison package?
  virtual std::vector<double> HandleEventRecord(const EventRecord *event) = 0;
  // virtual std::vector<double> GetObservablesValues() const = 0;
};
} // namespace rew
} // namespace genie

#endif