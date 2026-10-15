//____________________________________________________________________________
/*!

\class    genie::rew::RwgKineSpace

\brief    Abstract interface for algorithms that compute the kinematic
          variables and interaction channel ID used to classify an event for
          GENIE/Professor spline-based reweighting (see ObservableSplines).

\author   Qiyu Yan <yanqiyu17 \at mails.ucas.ac.cn>
          University of Chinese Academy of Sciences (UCAS) and University of Warwick

\created  July 2023

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

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
  virtual ChannelIDs ChannelID(const EventRecord &event) const = 0;

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
