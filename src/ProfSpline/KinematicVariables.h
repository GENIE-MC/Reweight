//____________________________________________________________________________
/*!

\class    genie::rew::ChannelIDs

\brief    A lightweight container identifying an interaction channel (e.g. by
          probe and target/hit-nucleon PDG codes) used to tag GENIE/Professor
          reweighting bins and splines.
          This header also defines genie::rew::KinematicVariables, a plain
          vector-of-doubles alias used to pass kinematic variable values
          around the ProfSpline machinery.

\author   Qiyu Yan <yanqiyu17 \at mails.ucas.ac.cn>
          University of Chinese Academy of Sciences (UCAS) and University of Warwick

\created  April 2024

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

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