//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightINukeParamsExtra

\brief   Helper class for FSI reweighting

\author  Gray Putnam <gputnam \at fnal.gov>
         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

\created Oct 8, 2025

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INTRANUKE_EXTRA_PARAMS_H_
#define _G_REWEIGHT_INTRANUKE_EXTRA_PARAMS_H_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINukeParams.h"
#include "RwCalculators/GReWeightINukeModelSwitchData.h"

namespace genie {
namespace rew   {

 class GReWeightINukeParamsExtra : public GReWeightINukeParams {

 public:

   GReWeightINukeParamsExtra();
  ~GReWeightINukeParamsExtra();

   class Fates;
   class MFP;

   //.........................................................................
   //
   // nested class: Fates
   //
   //.........................................................................

   class Fates : public GReWeightINukeParams::Fates {
   public :

     Fates(HadronType_t hadtype = kRwINukeUndefined);
    ~Fates();

     double OneSigmaErr(GSyst_t syst, double KE) const override;
     bool   IsG4ModelTransform(GSyst_t syst) const;
     bool   IsINCLModelTransform(GSyst_t syst) const;
     bool   IsModelTransform(GSyst_t syst) const;
     void   SetSystKERange(GSyst_t syst);
     bool   IsInSystKERange(double KE) const;

     virtual void Reset(void) override;
     virtual void SetTwkDial(GSyst_t s, double val) override;
     virtual double ScaleFactor(GSyst_t s, double KE ) const override;

     ModelSwitch_t        fModelSwitch;
     double fSystKELow; ///< Lower limit in kinetic energy range for this systematic
     double fSystKEHigh; ///< Upper limit in kinetic energy range for this systematic

   }; // Fates nested class


   //.........................................................................
   //
   // nested class: MFP
   //
   //.........................................................................

   class MFP : public GReWeightINukeParams::MFP {
   public :
     MFP( HadronType_t hadtype = kRwINukeUndefined );
     ~MFP();

     virtual void Reset(void) override;
     virtual double ScaleFactor( double KE ) const override;
     virtual void SetTwkDial( GSyst_t syst, double val ) override;
     void SetSystKERange( GSyst_t syst );
     bool IsInSystKERange( double KE ) const;

   private:
     double fSystKELow;
     double fSystKEHigh;

   }; // MFP nested class

 }; //GReWeightINukeParamsExtra

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_INTRANUKE_EXTRA_PARAMS_H_
