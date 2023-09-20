#ifndef _ObservablePMuEnuW_H_
#define _ObservablePMuEnuW_H_
#include "ProfSpline/ObservableI.h"

namespace genie {
namespace rew {
class ObservablePMuEnuW : public ObservableI {
public:
  virtual std::vector<double>
  GetKinematicVariables(const EventRecord &event) const override;

  virtual ~ObservablePMuEnuW() = default;
};
} // namespace rew
} // namespace genie

#endif