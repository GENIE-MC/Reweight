// Class to handle all observables as a whole
// This is a interface class

#ifndef _OBSERVABLE_I_
#define _OBSERVABLE_I_
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/EventGen/EventRecord.h"
#include <string>
#include <vector>
namespace genie {
namespace rew {
/// \uml{note Made this an interface to handle the possibliy of
/// different choice of different observables.
/// }
class ObservableI : public Algorithm {
public:
  // virtual void Configure(/*some way of configurating*/) = 0;

  ObservableI() = default;
  ObservableI(std::string name) : Algorithm(name) {}

  // calculate the value for each observable
  // maybe we can merge with binning lookup function below
  // but seperating them may benefit the idea of merging implementation
  // of corresponding ObservablePrediction in comparison package?
  virtual std::vector<double>
  GetKinematicVariables(const EventRecord &event) const = 0;

  virtual bool IsHandled(const EventRecord &event) const = 0;
  // IsCC is seperately exposed as it is crucial informaion 
  // of normalization, and I suppose that we will never want 
  // one ObservableI instance, and namely one GeneralReweightObs
  // instance to handle a mixture of CC and NC events
  //
  // So putting isCC a seperate flag to ensure that
  virtual bool IsCC() const = 0;
  // virtual std::vector<double> GetObservablesValues() const = 0;

  virtual ~ObservableI() = default;
};
} // namespace rew
} // namespace genie

#endif