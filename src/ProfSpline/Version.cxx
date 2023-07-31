#ifndef PROF_VERSION
#define PROF_VERSION 2.X
#endif

#define XSTR(s) STR(s)
#define STR(s) #s

#include <string>

namespace Professor {

  /// Version code for this build of Professor
  std::string version() {
    return XSTR(PROF_VERSION);
  }

}
