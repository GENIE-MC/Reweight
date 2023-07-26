#include "ProfSpline/ExtendedSpline.h"

namespace genie {
namespace rew {
ExtendedSpline::ExtendedSpline(const std::vector<double> &n_coefficients)
    : coefficients(n_coefficients) {}
ExtendedSpline::ExtendedSpline(const std::vector<double> &&n_coefficients)
    : coefficients(n_coefficients) {}
double ExtendedSpline::operator()(const std::vector<double> &) const {
  // TODO: we don't know what to do unless we understand how to interpret the
  // Prof2 spline

  // what we are supposed to do is to take the coefficients and the x values
  // to calculate the interpolated differential cross section
  return 0;
}
} // namespace rew
} // namespace genie