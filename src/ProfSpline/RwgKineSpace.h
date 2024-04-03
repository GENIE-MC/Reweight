// Class to handle all observables as a whole
// This is a interface class

#ifndef _OBSERVABLE_I_
#define _OBSERVABLE_I_
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/EventGen/EventRecord.h"
#include "ProfSpline/KinematicVariables.h"
#include <string>
#include <vector>
namespace genie {
namespace rew {
/// \uml{note Made this an interface to handle the possibliy of
/// different choice of different observables.
/// }
class RwgKineSpace : public Algorithm {
public:
  // calculate the value for each observable
  // maybe we can merge with binning lookup function below
  // but seperating them may benefit the idea of merging implementation
  // of corresponding ObservablePrediction in comparison package?
  virtual KinematicVariables
  CalcKinematicVariables(const EventRecord &event) const = 0;

  virtual bool IsHandled(const EventRecord &event) const = 0;
  // virtual std::vector<double> GetObservablesValues() const = 0;

  virtual void Configure(const Registry &config) override;
  virtual void Configure(string param_set) override;

  virtual ~RwgKineSpace() = default;

protected:
  virtual void LoadConfig(void) = 0;
  using Algorithm::Algorithm;
};
} // namespace rew
} // namespace genie

#endif