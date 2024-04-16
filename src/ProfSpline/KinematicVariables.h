#ifndef _KINEMATIC_VARIABLES_
#define _KINEMATIC_VARIABLES_
#include <algorithm>
#include <string>
#include <vector>

namespace genie {
namespace rew {

class ChannelIDs : public std::vector<int> {
public:
  using std::vector<int>::vector;
  ChannelIDs() = default;
  bool operator==(const ChannelIDs &o) const {
    return size() == o.size() && std::equal(begin(), end(), o.begin());
  }
  bool operator<(const ChannelIDs &o) const {
    return std::lexicographical_compare(begin(), end(), o.begin(), o.end());
  }
  bool operator>(const ChannelIDs &o) const { return o < *this; }
  bool operator<=(const ChannelIDs &o) const { return !(o < *this); }
  bool operator>=(const ChannelIDs &o) const { return !(*this < o); }
  bool operator!=(const ChannelIDs &o) const { return !(*this == o); }
  std::string toString(std::string delimiter = ",") const {
    std::string str = "";
    for (auto &i : *this) {
      str += std::to_string(i) + delimiter;
    }
    // remove the last delimiter
    str = str.substr(0, str.size() - delimiter.size());
    return str;
  }
  ChannelIDs(std::string str_in, std::string delimiter = ",") {
    std::string str = str_in;
    while (str.find(delimiter) != std::string::npos) {
      auto pos = str.find(delimiter);
      push_back(std::stoi(str.substr(0, pos)));
      str = str.substr(pos + 1);
    }
    push_back(std::stoi(str));
  }
};

using KinematicVariables = std::vector<double>;

// Just a wrapper class to hold the kinematic variables and the channel
// information, make it header only?
// class KinematicVariables {
// public:
//   KinematicVariables(std::vector<double> vars, ChannelIDs channel)
//       : vars(vars), channel(channel) {}

//   KinematicVariables() = default;
//   KinematicVariables(KinematicVariables &) = default;
//   KinematicVariables(KinematicVariables &&) = default;
//   KinematicVariables &operator=(KinematicVariables &) = default;
//   KinematicVariables &operator=(KinematicVariables &&) = default;
//   ~KinematicVariables() = default;

//   std::vector<double> &GetVars() { return vars; }
//   ChannelIDs &GetChannel() { return channel; }

//   const std::vector<double> &GetVars() const { return vars; }
//   const ChannelIDs &GetChannel() const { return channel; }

//   bool MatchChannel(const ChannelIDs &o_channel) const {
//     return channel.size() == o_channel.size() &&
//            std::equal(channel.begin(), channel.end(), o_channel.begin());
//   }

// private:
//   std::vector<double> vars;
//   ChannelIDs channel;
// };

} // namespace rew
} // namespace genie

#endif