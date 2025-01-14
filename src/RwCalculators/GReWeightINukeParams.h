//____________________________________________________________________________
/*!

\class   genie::rew::GReWeightINukeParams

\brief   Helper class for cross section model reweighting

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
         University of Liverpool

         Jim Dobson <J.Dobson07 \at imperial.ac.uk>
         Imperial College London

\created Sep 10, 2009

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_INTRANUKE_PARAMS_H_
#define _G_REWEIGHT_INTRANUKE_PARAMS_H_

#include <map>

// GENIE/Generator includes
#include "Framework/ParticleData/PDGUtils.h"

// GENIE/Reweight includes
#include "RwFramework/GSyst.h"

class TLorentzVector;

namespace genie {
namespace rew   {

 class GReWeightINukeParams {

 public:

   typedef enum EHadronType {
      kRwINukeUndefined = 0,
      kRwINukePion,
      kRwINukeNucl,
   } HadronType_t;

   static HadronType_t HadronTypeFromPdg(int pdgc) {
     if(pdg::IsPion   (pdgc)) return kRwINukePion;
     if(pdg::IsNucleon(pdgc)) return kRwINukeNucl;
     return kRwINukeUndefined;
   }

   GReWeightINukeParams();
  ~GReWeightINukeParams();

   class Fates;
   class MFP;

   Fates * FateParams         (int pdgc) const;         ///<
   MFP *   MeanFreePathParams (int pdgc) const;         ///<
   void    Reset              (void);                   ///<
   void    Reconfigure        (void);                   ///<
   double  ChisqPenalty       (void) const;             ///<
   void    SetTwkDial         (GSyst_t s, double val);  ///<
   void    SetTargetA         (int target_A); ///< Set the mass number of the hit nucleus

   //.........................................................................
   //
   // nested class: Fates
   //
   //.........................................................................

   class Fates {
   public :
     Fates(HadronType_t hadtype = kRwINukeUndefined);
    ~Fates();

     double ScaleFactor   (GSyst_t s, const TLorentzVector & p4) const; ///< see next
     double ScaleFactor   (GSyst_t s, double KE=-1.) const;             ///< fate fraction scale factor = 1 + twk_dial * fractional_err
     bool   IsIncluded    (GSyst_t s) const;                            ///< is included?
     bool   IsCushionTerm (GSyst_t s) const;                            ///< is it a cushion term?
     bool   IsTweaked     (GSyst_t s) const;                            ///< is included & tweaked to non-def value?
     bool   IsTweaked     (void) const;                                 ///< is any param tweaked
     void   Reset         (void);                                       ///<
     void   Reconfigure   (void);                                       ///<
     double ChisqPenalty  (void) const;                                 ///<
     void   SetTwkDial    (GSyst_t s, double val);                      ///<
     void   SetTargetA    (int target_A); ///< Set the mass number of the hit nucleus

   private:

     bool   IsHandled       (GSyst_t s) const;
     void   AddCushionTerms (void);
     double ActualTwkDial   (GSyst_t s, double KE=-1.) const;  ///< actual tweaking dial for input systematic at input kinetic energy

     HadronType_t         fHadType;           ///<
     std::map<GSyst_t, double> fSystValuesUser;    ///< List of systematics included & values set by the user
     mutable std:: map<GSyst_t, double> fSystValuesActual;  ///< List of systematics included & values actually used (user values limited to physical range)
     std::map<GSyst_t, bool>   fIsCushion;         ///< cushion term flag
     int fTargetA; ///< Mass number of the hit nucleus (needed for pion fates)

   }; // Fates nested class


   //.........................................................................
   //
   // nested class: MFP
   //
   //.........................................................................

   class MFP {
   public :
     MFP(HadronType_t hadtype = kRwINukeUndefined);
    ~MFP();

     double ScaleFactor   (void) const;  ///< mean free path scale factor = 1 + twk_dial * fractional_err
     double TwkDial       (void) const;  ///< current value of mfp tweak dial
     bool   IsIncluded    (void) const;  ///<
     bool   IsTweaked     (void) const;  ///<
     double ChisqPenalty  (void) const;  ///<
     void   Reset         (void);        ///<
     void   SetTwkDial    (double val);  ///<

   private:
     HadronType_t fHadType;     ///<
     GSyst_t      fSyst;        ///<
     double       fTwkDial;     ///<
     bool         fIsIncluded;  ///<

   }; // MFP nested class


 private:

    Fates * fParmPionFates;
    Fates * fParmNuclFates;
    MFP *   fParmPionMFP;
    MFP *   fParmNuclMFP;

 }; //GReWeightINukeParams

}      // rew   namespace
}      // genie namespace

#endif // _G_REWEIGHT_INTRANUKE_PARAMS_H_
