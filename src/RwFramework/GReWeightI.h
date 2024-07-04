//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightI

\brief    GENIE event reweighting engine ABC 

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_ABC_H_
#define _G_REWEIGHT_ABC_H_

// GENIE/Generator includes
#include "Framework/Interaction/ScatteringType.h"
#include "Framework/EventGen/EventRecord.h"

// GENIE/Reweight includes
#include "RwFramework/GSyst.h"
#include "RwFramework/GSystSet.h"

namespace genie {

class EventRecord;

namespace rew   {

 class GReWeightI 
 {
 public:
  virtual ~GReWeightI() 
  { 

  } 

  //
  // define the GReWeightI interface
  //
  
  //! does the current weight calculator handle this type of event?
  virtual bool AppliesTo (const EventRecord & event) const = 0; 

  //! does the current weight calculator handle the input nuisance param?
  virtual bool IsHandled (GSyst_t syst) const = 0;

  //! update the value for the specified nuisance param
  virtual void SetSystematic (GSyst_t syst, double val) = 0; 

  //!  set all nuisance parameters to default values
  virtual void Reset (void) = 0;            

  //! propagate updated nuisance parameter values to actual MC, etc
  virtual void Reconfigure (void) = 0;            
  
  //! calculate a weight for the input event using the current nuisance param values
  virtual double CalcWeight (const genie::EventRecord & event) = 0;
  
  //! Should we calculate the old weight ourselves, or use the one from the input tree? Default on.
  virtual void UseOldWeightFromFile(bool) = 0;
  
  //! If using the weight from the file, how many times should we check by calculating it ourself?
  virtual void SetNWeightChecks(int) = 0;

 protected:

   GReWeightI() 
   { 
   
   }
 };

} // rew   namespace
} // genie namespace

#endif

