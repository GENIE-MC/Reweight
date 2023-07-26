// Actual storage for bin information
// Should be initialized by ObservableSplines
//

#ifndef _EXTENDED_SPLINE_
#define _EXTENDED_SPLINE_
#include <vector>

namespace genie {
namespace rew {
class ExtendedSpline {
public:
  ExtendedSpline(const std::vector<double> &coefficients);
  ExtendedSpline(const std::vector<double> &&coefficients);
  double operator()(const std::vector<double> &x) const;

private:
  std::vector<double> coefficients; // to be read from Professor spline file
};
} // namespace rew
} // namespace genie

#endif