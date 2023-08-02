// Class to handle all observables as a whole
// This is a interface class

#ifndef _OBSERVABLE_I_
#define _OBSERVABLE_I_
#include "Framework/EventGen/EventRecord.h"
#include <string>
#include <vector>
namespace genie {
namespace rew {
/// \uml{note Made this an interface to handle the possibliy of
/// different choice of different observables.
///
/// TODO: consider how to initialize the Observable class
/// TODO: how to handle different set of observables in different
/// channels
/// IDEA: Let GReWeightProfessor hold multiple ObservableSplines?
/// TODO: Any implementation of ObservableI would do have corresponding
/// ObservablePrediction in comparison package, doing the same thing
/// Should we consider merging them together?
/// }
class ObservableI {
public:
  // virtual void Configure(/*some way of configurating*/) = 0;

  // calculate the value for each observable
  // maybe we can merge with binning lookup function below
  // but seperating them may benefit the idea of merging implementation
  // of corresponding ObservablePrediction in comparison package?
  virtual std::vector<double> GetKinematicVariables(const EventRecord &event) = 0;
  // virtual std::vector<double> GetObservablesValues() const = 0;

  virtual ~ObservableI() = default;
};
} // namespace rew
} // namespace genie

#endif