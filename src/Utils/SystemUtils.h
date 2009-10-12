//____________________________________________________________________________
/*!

\namespace  genie::utils::system

\brief      System utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            STFC, Rutherford Appleton Laboratory

\created    Oct 08, 2009

\cpright    Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SYST_UTILS_H_
#define _SYST_UTILS_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace genie {
namespace utils {

namespace system
{

  vector<string> GetAllFilesInPath  (string path, string extension="");

} // system namespace
} // utils  namespace
} // genie  namespace

#endif // _SYST_UTILS_H_