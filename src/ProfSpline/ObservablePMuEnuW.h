#ifndef _ObservablePMuEnuW_H_
#define _ObservablePMuEnuW_H_
#include "ProfSpline/RwgKineSpace.h"

namespace genie {
namespace rew {
class ObservablePMuEnuW : public RwgKineSpace {
public:
  virtual std::vector<double>
  KinematicVariables(const EventRecord &event) const override;

  virtual ~ObservablePMuEnuW() = default;
};
} // namespace rew
} // namespace genie

#endif