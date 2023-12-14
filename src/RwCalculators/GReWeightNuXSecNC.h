//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecNC

\brief    Tweak NC level in a gross way (scale by a given factor).
          More detailed rew. functions also exist for specific NC sources.
          RES NC cross section rew. is also handled by GReWeightNuXSecNCRES.
          DIS NC cross section rew. is also handled by GReWeightNuXSecDIS.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  May 25, 2010

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_NC_H_
#define _G_REWEIGHT_NU_XSEC_NC_H_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

namespace genie {
namespace rew   {

 class GReWeightNuXSecNC : public GReWeightModel
 {
 public:
   GReWeightNuXSecNC();
  ~GReWeightNuXSecNC();

   // implement the GReWeightI interface
   bool   AppliesTo      (const EventRecord & event) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void RewNue       (bool   tf)  { fRewNue     = tf; }
   void RewNuebar    (bool   tf)  { fRewNuebar  = tf; }
   void RewNumu      (bool   tf)  { fRewNumu    = tf; }
   void RewNumubar   (bool   tf)  { fRewNumubar = tf; }
   void RewQE        (bool   tf)  { fRewQE      = tf; }
   void RewRES       (bool   tf)  { fRewRES     = tf; }
   void RewDIS       (bool   tf)  { fRewDIS     = tf; }

 private:

   void   Init (void);

   bool   fRewNue;         ///< reweight nu_e?
   bool   fRewNuebar;      ///< reweight nu_e_bar?
   bool   fRewNumu;        ///< reweight nu_mu?
   bool   fRewNumubar;     ///< reweight nu_mu_bar?
   bool   fRewQE;          ///< reweight NC elastic?
   bool   fRewRES;         ///< reweight NC RES?
   bool   fRewDIS;         ///< reweight NC DIS?
   double fNCTwkDial;      ///< tweaking dial

   std::string fManualModelName; ///< If using a tweaked model that isn't the same as default, name
   std::string fManualModelType; ///< If using a tweaked model that isn't the same as default, type
 };

} // rew   namespace
} // genie namespace

#endif
