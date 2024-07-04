//____________________________________________________________________________
/*!

\class    genie::rew::GSystSet

\brief    Set of systematics to be considered by the reweighting package.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_SET_OF_SYSTEMATICS_H_
#define _G_SET_OF_SYSTEMATICS_H_

#include <string>
#include <map>
#include <vector>

// GENIE/Reweight includes
#include "RwFramework/GSyst.h"

using std::string;
using std::map;
using std::vector;

namespace genie {
namespace rew   {

class GSystInfo;

class GSystSet {

public:  
  GSystSet();
  GSystSet(const GSystSet & syst_set);
 ~GSystSet();

  void Init    (GSyst_t syst, double init=0., double min=-1., double max=+1., double step=0.05);
  void Remove  (GSyst_t syst);
  void Set     (GSyst_t syst, double current_value);
  int  Size    (void) const;
  bool Added   (GSyst_t syst) const;
  void Print   (void);
  void Copy    (const GSystSet & syst_set);

  const GSystInfo * Info(GSyst_t syst) const;

  vector<genie::rew::GSyst_t> AllIncluded (void);

private:
  
  map<GSyst_t, GSystInfo *>  fSystematics;  
};

class GSystInfo {

public:
  GSystInfo() : 
     CurValue(0), InitValue(0), MinValue(0), MaxValue(0), Step(0) 
  { 

  }
  GSystInfo(double init, double min, double max, double step) : 
     CurValue(init), InitValue(init), MinValue(min), MaxValue(max), Step(step) 
  {
 
  }
 ~GSystInfo() 
  { 

  }

  double CurValue;
  double InitValue;
  double MinValue;
  double MaxValue;
  double Step;
};

} // rew   namespace
} // genie namespace

#endif 

