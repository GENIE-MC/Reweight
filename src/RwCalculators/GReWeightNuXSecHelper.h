//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightNuXSecHelper

\brief   Helper class for cross section model reweighting

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
         University of Liverpool

\created October 22, 2005

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NEUTRINO_CROSS_SECTION_HELPER_H_
#define _G_REWEIGHT_NEUTRINO_CROSS_SECTION_HELPER_H_

#include <map>

// GENIE/Generator includes
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Interaction/ScatteringType.h"
#include "Framework/EventGen/GEVGPool.h"

namespace genie {

class EventRecord;

namespace rew {

class GReWeightNuXSecHelper {

public :
  GReWeightNuXSecHelper();
 ~GReWeightNuXSecHelper();

  void   HandleInitState  (const InitialState & init_state);
  void   DiffCrossSecType (ScatteringType_t sct, KinePhaseSpace_t kps);
  double NewWeight        (const EventRecord & event, bool shape_only = false);

private:

   void Initialize();

   GEVGPool                                fGPool;             ///<
   std::map<ScatteringType_t, KinePhaseSpace_t> fCrossSecModelPhSp; ///<
};

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_NEUTRINO_CROSS_SECTION_H_
