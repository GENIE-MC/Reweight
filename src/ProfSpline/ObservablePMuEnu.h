#ifndef _ObservablePMuEnu_H_
#define _ObservablePMuEnu_H_
#include "ProfSpline/ObservableI.h"

namespace genie {
namespace rew {
class ObservablePMuEnu : public RwgKineSpace {
public:
  virtual std::vector<double>
  KinematicVariables(const EventRecord &event) const override;

  virtual ~ObservablePMuEnu() = default;
};
} // namespace rew
} // namespace genie

#endif