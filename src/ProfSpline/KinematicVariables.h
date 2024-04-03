#ifndef _KINEMATIC_VARIABLES_
#define _KINEMATIC_VARIABLES_
#include <algorithm>
#include <vector>

namespace genie {
namespace rew {

// Just a wrapper class to hold the kinematic variables and the channel
// information, make it header only?
class KinematicVariables {
public:
  KinematicVariables(std::vector<double> vars, std::vector<int> channel)
      : vars(vars), channel(channel) {}

  KinematicVariables() = default;
  KinematicVariables(KinematicVariables &) = default;
  KinematicVariables(KinematicVariables &&) = default;
  KinematicVariables &operator=(KinematicVariables &) = default;
  KinematicVariables &operator=(KinematicVariables &&) = default;
  ~KinematicVariables() = default;

  std::vector<double> &GetVars() { return vars; }
  std::vector<int> &GetChannel() { return channel; }

  const std::vector<double> &GetVars() const { return vars; }
  const std::vector<int> &GetChannel() const { return channel; }

  bool MatchChannel(const std::vector<int> &o_channel) const {
    return channel.size() == o_channel.size() &&
           std::equal(channel.begin(), channel.end(), o_channel.begin());
  }

private:
  std::vector<double> vars;
  std::vector<int> channel;
};

} // namespace rew
} // namespace genie

#endif