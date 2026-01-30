//____________________________________________________________________________
/*!

\class    genie::rew::GSyst_t

\brief    An enumeration of systematic parameters

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_SYSTEMATIC_PARAM_H_
#define _G_SYSTEMATIC_PARAM_H_

#include <map>
#include <string>

// GENIE/Generator includes
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Interaction/InteractionType.h"
#include "Physics/HadronTransport/INukeHadroFates2018.h"

using std::string;

namespace genie {
namespace rew   {

typedef enum EGSyst {

  kNullSystematic = 0,

  //
  // Neutrino cross section systematics
  //

  // NCEL tweaking parameters:
  kXSecTwkDial_MaNCEL,            ///< tweak Ma NCEL, affects dsigma(NCEL)/dQ2 both in shape and normalization
  kXSecTwkDial_EtaNCEL,           ///< tweak NCEL strange axial form factor eta, affects dsigma(NCEL)/dQ2 both in shape and normalization
  // CCQE tweaking parameters (also see at end of enum for z-expansion FF knobs):
  kXSecTwkDial_NormCCQE,          ///< tweak CCQE normalization (energy independent)
  kXSecTwkDial_NormCCQEenu,       ///< tweak CCQE normalization (maintains dependence on neutrino energy)
  kXSecTwkDial_MaCCQEshape,       ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCQE,            ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_VecFFCCQEshape,    ///< tweak elastic nucleon form factors (BBA/default -> dipole) - shape only effect of dsigma(CCQE)/dQ2
  // Resonance neutrino-production tweaking parameters:
  kXSecTwkDial_NormCCRES,         ///< tweak CCRES normalization
  kXSecTwkDial_MaCCRESshape,      ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MvCCRESshape,      ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCRES,           ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MvCCRES,           ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_NormNCRES,         ///< tweak NCRES normalization
  kXSecTwkDial_MaNCRESshape,      ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MvNCRESshape,      ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaNCRES,           ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MvNCRES,           ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
  // Coherent pion production tweaking parameters:
  kXSecTwkDial_MaCOHpi,           ///< tweak Ma for COH pion production
  kXSecTwkDial_R0COHpi,           ///< tweak R0 for COH pion production
  // Non-resonance background tweaking parameters:
  kXSecTwkDial_RvpCC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+p CC
  kXSecTwkDial_RvpCC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+p CC
  kXSecTwkDial_RvpNC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+p NC
  kXSecTwkDial_RvpNC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+p NC
  kXSecTwkDial_RvnCC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+n CC
  kXSecTwkDial_RvnCC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+n CC
  kXSecTwkDial_RvnNC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+n NC
  kXSecTwkDial_RvnNC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+n NC
  kXSecTwkDial_RvbarpCC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p CC
  kXSecTwkDial_RvbarpCC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p CC
  kXSecTwkDial_RvbarpNC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p NC
  kXSecTwkDial_RvbarpNC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p NC
  kXSecTwkDial_RvbarnCC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n CC
  kXSecTwkDial_RvbarnCC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n CC
  kXSecTwkDial_RvbarnNC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n NC
  kXSecTwkDial_RvbarnNC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n NC
  // DIS tweaking parameters - applied for DIS events with (Q2>Q2o, W>Wo), typically Q2o=1GeV^2, Wo=1.7-2.0GeV
  kXSecTwkDial_AhtBY,             ///< tweak the Bodek-Yang model parameter A_{ht} - incl. both shape and normalization effect
  kXSecTwkDial_BhtBY,             ///< tweak the Bodek-Yang model parameter B_{ht} - incl. both shape and normalization effect
  kXSecTwkDial_CV1uBY,            ///< tweak the Bodek-Yang model parameter CV1u - incl. both shape and normalization effect
  kXSecTwkDial_CV2uBY,            ///< tweak the Bodek-Yang model parameter CV2u - incl. both shape and normalization effect
  kXSecTwkDial_AhtBYshape,        ///< tweak the Bodek-Yang model parameter A_{ht} - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_BhtBYshape,        ///< tweak the Bodek-Yang model parameter B_{ht} - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_CV1uBYshape,       ///< tweak the Bodek-Yang model parameter CV1u - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_CV2uBYshape,       ///< tweak the Bodek-Yang model parameter CV2u - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_NormDISCC,         ///< tweak the inclusive DIS CC normalization
  kXSecTwkDial_RnubarnuCC,        ///< tweak the ratio of \sigma(\bar\nu CC) / \sigma(\nu CC)
  kXSecTwkDial_DISNuclMod,        ///< tweak DIS nuclear modification (shadowing, anti-shadowing, EMC)
  //
  kXSecTwkDial_NC,                ///<


  //
  // Hadronization (free nucleon target)
  //

  kHadrAGKYTwkDial_xF1pi,         ///< tweak xF distribution for low multiplicity (N + pi) DIS f/s produced by AGKY
  kHadrAGKYTwkDial_pT1pi,         ///< tweak pT distribution for low multiplicity (N + pi) DIS f/s produced by AGKY

  //
  // Medium-effects to hadronization
  //

  kHadrNuclTwkDial_FormZone,         ///< tweak formation zone

  //
  // Intranuclear rescattering systematics.
  // There are 2 sets of parameters:
  // - parameters that control the total rescattering probability, P(total)
  // - parameters that control the fraction of each process (`fate'), given a total rescat. prob., P(fate|total)
  // These parameters are considered separately for pions and nucleons.
  //

  kINukeTwkDial_MFP_pi,      ///< tweak mean free path for pions
  kINukeTwkDial_MFP_N,       ///< tweak mean free path for nucleons
  kINukeTwkDial_MFPLoE_N,       ///< tweak mean free path for nucleons, 0 <= KE < 150 MeV
  kINukeTwkDial_MFPM1E_N,       ///< tweak mean free path for nucleons, 150 <= KE < 300 MeV
  kINukeTwkDial_MFPM2E_N,       ///< tweak mean free path for nucleons, 300 <= KE < 600 MeV
  kINukeTwkDial_MFPHiE_N,       ///< tweak mean free path for nucleons, 600 <= KE MeV
  kINukeTwkDial_FrCEx_pi,    ///< tweak charge exchange probability for pions, for given total rescattering probability
  // sd - hA no longer has elastic fate
  kINVALID_INukeTwkDial_FrElas_pi,   ///< tweak elastic         probability for pions, for given total rescattering probability
  kINukeTwkDial_FrInel_pi,   ///< tweak inelastic       probability for pions, for given total rescattering probability
  kINukeTwkDial_FrAbs_pi,    ///< tweak absorption      probability for pions, for given total rescattering probability
  kINukeTwkDial_FrPiProd_pi, ///< tweak pion production probability for pions, for given total rescattering probability

  kINukeTwkDial_G4_N,     ///< tweak intranuclear scattering to G4 values
  kINukeTwkDial_INCL_N,     ///< tweak intranuclear scattering to INCL values
  kINukeTwkDial_G4LoE_N,     ///< tweak intranuclear scattering to G4 values, 0 <= KE < 150 MeV
  kINukeTwkDial_INCLLoE_N,     ///< tweak intranuclear scattering to INCL values, 0 <= KE < 150 MeV
  kINukeTwkDial_G4M1E_N,     ///< tweak intranuclear scattering to G4 values, 150 <= KE < 300 MeV
  kINukeTwkDial_INCLM1E_N,     ///< tweak intranuclear scattering to INCL values, 300 <= KE < 300 MeV
  kINukeTwkDial_G4M2E_N,     ///< tweak intranuclear scattering to G4 values, 300 <= KE < 600 MeV
  kINukeTwkDial_INCLM2E_N,     ///< tweak intranuclear scattering to INCL values, 300 <= KE < 600 MeV
  kINukeTwkDial_G4HiE_N,     ///< tweak intranuclear scattering to G4 values, 600 <= KE
  kINukeTwkDial_INCLHiE_N,     ///< tweak intranuclear scattering to INCL values, 600 <= KE

  kINukeTwkDial_FrCEx_N,     ///< tweak charge exchange probability for nucleons, for given total rescattering probability
  kINVALID_INukeTwkDial_FrElas_N,    ///< tweak elastic         probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrInel_N,    ///< tweak inelastic       probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrAbs_N,     ///< tweak absorption      probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrPiProd_N,  ///< tweak pion production probability for nucleons, for given total rescattering probability

  //
  // Nuclear model
  //

  kSystNucl_CCQEPauliSupViaKF,   ///<
  kSystNucl_CCQEMomDistroFGtoSF, ///<

  //
  // Resonance decays
  //

  kRDcyTwkDial_BR1gamma,        ///< tweak Resonance -> X + gamma branching ratio, eg Delta+(1232) -> p gamma
  kRDcyTwkDial_BR1eta,          ///< tweak Resonance -> X + eta   branching ratio, eg N+(1440) -> p eta
  kRDcyTwkDial_Theta_Delta2Npi,  ///< distort pi angular distribution in Delta -> N + pi

  //
  // Alternative approach to CCQE form factors (z-expansion)
  //

  kXSecTwkDial_ZNormCCQE,         ///< tweak Z-expansion CCQE normalization (energy independent)
  kXSecTwkDial_ZExpA1CCQE,        ///< tweak Z-expansion coefficient 1, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_ZExpA2CCQE,        ///< tweak Z-expansion coefficient 2, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_ZExpA3CCQE,        ///< tweak Z-expansion coefficient 3, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_ZExpA4CCQE,        ///< tweak Z-expansion coefficient 4, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_AxFFCCQEshape,     ///< tweak axial nucleon form factors (dipole -> z-expansion) - shape only effect of dsigma(CCQE)/dQ2

  //
  //   Alternative approach to CCQE form factors (RunningMA)
  //

  kXSecTwkDial_E0CCQEshape,       ///< tweak E0 CCQE RunningMA, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_E0CCQE,            ///< tweak E0 CCQE RunningMA, affects dsigma(CCQE)/dQ2 both in shape and normalization

  //
  // Empirical MEC dials
  //
  kXSecTwkDial_EmpMEC_Mq2d,
  kXSecTwkDial_EmpMEC_Mass,
  kXSecTwkDial_EmpMEC_Width,
  kXSecTwkDial_EmpMEC_FracPN_NC,
  kXSecTwkDial_EmpMEC_FracPN_CC,
  kXSecTwkDial_EmpMEC_FracCCQE,
  kXSecTwkDial_EmpMEC_FracNCQE,
  kXSecTwkDial_EmpMEC_FracPN_EM,
  kXSecTwkDial_EmpMEC_FracEMQE,

  //
  // General MEC dials
  //

  // MEC cross section normalization
  kXSecTwkDial_NormCCMEC,
  kXSecTwkDial_NormNCMEC,
  kXSecTwkDial_NormEMMEC,

  // MEC nucleon cluster decay angular distribution
  kXSecTwkDial_DecayAngMEC,

  // Fraction of CCMEC initial nucleon clusters that are p+n
  kXSecTwkDial_FracPN_CCMEC,

  // Fraction of CCMEC events involving an internal delta line
  // (currently only used by Valencia MEC)
  kXSecTwkDial_FracDelta_CCMEC,

  // Shape of CCMEC differential cross section (interpolates
  // between models)
  kXSecTwkDial_XSecShape_CCMEC,

  // Interpolates between the default CCQE model and the same
  // one with RPA off (only gives non-unit weights for Nieves CCQE)
  kXSecTwkDial_RPA_CCQE,

  /// Distort photon angular distribution in Delta -> N + photon
  kRDcyTwkDial_Theta_Delta2NRad,

  // Tweak the value of the EM potential used when computing the Coulomb
  // correction factor in the Nieves CCQE model
  kXSecTwkDial_CoulombCCQE,

  // Scale the normalization of CC coherent pion production
  kXSecTwkDial_NormCCCOHpi,

  // Scale the normalization of NC coherent pion production
  kXSecTwkDial_NormNCCOHpi,

  //
  // Alternative approach to CCQE form factors (z-expansion) vector form factor
  //
  kXSecTwkDial_ZExpELFF,
  kXSecTwkDial_ZExpELFF_AP1,
  kXSecTwkDial_ZExpELFF_AP2,
  kXSecTwkDial_ZExpELFF_AP3,
  kXSecTwkDial_ZExpELFF_AP4,
  kXSecTwkDial_ZExpELFF_AN1,
  kXSecTwkDial_ZExpELFF_AN2,
  kXSecTwkDial_ZExpELFF_AN3,
  kXSecTwkDial_ZExpELFF_AN4,
  kXSecTwkDial_ZExpELFF_BP1,
  kXSecTwkDial_ZExpELFF_BP2,
  kXSecTwkDial_ZExpELFF_BP3,
  kXSecTwkDial_ZExpELFF_BP4,
  kXSecTwkDial_ZExpELFF_BN1,
  kXSecTwkDial_ZExpELFF_BN2,
  kXSecTwkDial_ZExpELFF_BN3,
  kXSecTwkDial_ZExpELFF_BN4,

  //
  // Misc
  //

  kNTwkDials, /// < Not a real dial, just keep as last entry for looping purposes

  kProfRew    /// Systematics to be handled by the Professor reweighting framework
} GSyst_t;


class GSyst {

public:
 //......................................................................................
 static string AsString(GSyst_t syst) {
   // Convert the GSyst_t value into its corresponding string by looking it
   // up in the map.
   std::map<GSyst_t, string>::const_iterator it = fGSystToStringMap.find( syst );
   if ( it != fGSystToStringMap.cend() ) return it->second;
   // If a match could not be found, then return "-".
   else return "-";
 }
 //......................................................................................
 static GSyst_t FromString(string syst_name) {
   // This search could be simplified a bit in C++11 (e.g., using std::find_if
   // and a lambda). Perhaps we could tweak it for GENIE 4. - S. Gardiner
   std::map<GSyst_t, string>::const_iterator it = fGSystToStringMap.cbegin();
   std::map<GSyst_t, string>::const_iterator end = fGSystToStringMap.cend();
   while ( it != end ) {
     // If the value stored in the map (i.e., the string containing the
     // systematic tweak knob name) matches the input string, return the
     // corresponding key (GSyst_t value)
     if ( it->second == syst_name ) return it->first;
     ++it;
   }

   // A matching string could not be found in the map, so just return the null
   // enum value
   return kNullSystematic;
 }
 //......................................................................................
 static bool IsINukePionFateSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_FrCEx_pi   ) :
       //     case ( kINukeTwkDial_FrElas_pi  ) :
     case ( kINukeTwkDial_FrInel_pi  ) :
     case ( kINukeTwkDial_FrAbs_pi   ) :
     case ( kINukeTwkDial_FrPiProd_pi) :
        return true;
        break;
     default:
        return false;
        break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeNuclFateSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_FrCEx_N   ) :
       //     case ( kINukeTwkDial_FrElas_N  ) :
     case ( kINukeTwkDial_FrInel_N  ) :
     case ( kINukeTwkDial_FrAbs_N   ) :
     case ( kINukeTwkDial_FrPiProd_N) :
     case ( kINukeTwkDial_G4_N) :
     case ( kINukeTwkDial_INCL_N) :
     case ( kINukeTwkDial_G4LoE_N) :
     case ( kINukeTwkDial_INCLLoE_N) :
     case ( kINukeTwkDial_G4M1E_N) :
     case ( kINukeTwkDial_INCLM1E_N) :
     case ( kINukeTwkDial_G4M2E_N) :
     case ( kINukeTwkDial_INCLM2E_N) :
     case ( kINukeTwkDial_G4HiE_N) :
     case ( kINukeTwkDial_INCLHiE_N) :
        return true;
        break;
     default:
        return false;
        break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeFateSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_FrCEx_pi    ) :
       //     case ( kINukeTwkDial_FrElas_pi   ) :
     case ( kINukeTwkDial_FrInel_pi   ) :
     case ( kINukeTwkDial_FrAbs_pi    ) :
     case ( kINukeTwkDial_FrPiProd_pi ) :
     case ( kINukeTwkDial_FrCEx_N     ) :
       //     case ( kINukeTwkDial_FrElas_N    ) :
     case ( kINukeTwkDial_FrInel_N    ) :
     case ( kINukeTwkDial_FrAbs_N     ) :
     case ( kINukeTwkDial_FrPiProd_N  ) :
       return true;
       break;

     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukePionMeanFreePathSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_MFP_pi ) :
       return true;
       break;

     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeNuclMeanFreePathSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_MFP_N  ) :
     case ( kINukeTwkDial_MFPLoE_N  ) :
     case ( kINukeTwkDial_MFPM1E_N  ) :
     case ( kINukeTwkDial_MFPM2E_N  ) :
     case ( kINukeTwkDial_MFPHiE_N  ) :
       return true;
       break;

     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeMeanFreePathSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_MFP_pi ) :
     case ( kINukeTwkDial_MFP_N  ) :
     case ( kINukeTwkDial_MFPLoE_N  ) :
     case ( kINukeTwkDial_MFPM1E_N  ) :
     case ( kINukeTwkDial_MFPM2E_N  ) :
     case ( kINukeTwkDial_MFPHiE_N  ) :
       return true;
       break;

     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static GSyst_t NextPionFateSystematic(int i)
 {
    if(i==0) return kINukeTwkDial_FrCEx_pi;
    //    if(i==1) return kINukeTwkDial_FrElas_pi;
    if(i==1) return kINukeTwkDial_FrInel_pi;
    if(i==2) return kINukeTwkDial_FrAbs_pi;
    if(i==3) return kINukeTwkDial_FrPiProd_pi;

    return kNullSystematic;
 }
 //......................................................................................
 static GSyst_t NextNuclFateSystematic(int i)
 {
    if(i==0) return kINukeTwkDial_FrCEx_N;
    //    if(i==1) return kINukeTwkDial_FrElas_N;
    if(i==1) return kINukeTwkDial_FrInel_N;
    if(i==2) return kINukeTwkDial_FrAbs_N;
    if(i==3) return kINukeTwkDial_FrPiProd_N;

    return kNullSystematic;
 }
 //......................................................................................
 static GSyst_t INukeFate2GSyst(INukeFateHA_t fate, int pdgc)
 {
  // get the corresponding GSyst_t systematic parameter enumeration from the
  // input intranuke fate enumeration and PDG code
  //
  if(pdg::IsPion(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kNullSystematic;            break;
      case kIHAFtCEx       : return kINukeTwkDial_FrCEx_pi;     break;
	//      case kIHAFtElas      : return kINukeTwkDial_FrElas_pi;    break;
      case kIHAFtInelas    : return kINukeTwkDial_FrInel_pi;    break;
      case kIHAFtAbs       : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtPiProd    : return kINukeTwkDial_FrPiProd_pi;  break;
      default              : return kNullSystematic;            break;
     }
  } else
  if(pdg::IsNucleon(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kNullSystematic;           break;
      case kIHAFtCEx       : return kINukeTwkDial_FrCEx_N;     break;
	//      case kIHAFtElas      : return kINukeTwkDial_FrElas_N;    break;
      case kIHAFtInelas    : return kINukeTwkDial_FrInel_N;    break;
      case kIHAFtAbs       : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtPiProd    : return kINukeTwkDial_FrPiProd_N;  break;
      default              : return kNullSystematic;           break;
     }
  }
  return kNullSystematic;
 }
 //......................................................................................
 static GSyst_t RBkg(InteractionType_t itype, int probe, int hitnuc, int npi)
 {
   bool is_v    = pdg::IsNeutrino     (probe);
   bool is_vbar = pdg::IsAntiNeutrino (probe);
   bool is_p    = pdg::IsProton       (hitnuc);
   bool is_n    = pdg::IsNeutron      (hitnuc);

   // CC
   bool is_cc = (itype == kIntWeakCC);
   if(is_cc) {
     if(is_v && is_p) {
       if(npi==1) return kXSecTwkDial_RvpCC1pi;
       if(npi==2) return kXSecTwkDial_RvpCC2pi;
     }
     if(is_v && is_n) {
       if(npi==1) return kXSecTwkDial_RvnCC1pi;
       if(npi==2) return kXSecTwkDial_RvnCC2pi;
     }
     if(is_vbar && is_p) {
       if(npi==1) return kXSecTwkDial_RvbarpCC1pi;
       if(npi==2) return kXSecTwkDial_RvbarpCC2pi;
     }
     if(is_vbar && is_n) {
       if(npi==1) return kXSecTwkDial_RvbarnCC1pi;
       if(npi==2) return kXSecTwkDial_RvbarnCC2pi;
     }
   }//cc

   // NC
   bool is_nc = (itype == kIntWeakNC);
   if(is_nc) {
     if(is_v && is_p) {
       if(npi==1) return kXSecTwkDial_RvpNC1pi;
       if(npi==2) return kXSecTwkDial_RvpNC2pi;
     }
     if(is_v && is_n) {
       if(npi==1) return kXSecTwkDial_RvnNC1pi;
       if(npi==2) return kXSecTwkDial_RvnNC2pi;
     }
     if(is_vbar && is_p) {
       if(npi==1) return kXSecTwkDial_RvbarpNC1pi;
       if(npi==2) return kXSecTwkDial_RvbarpNC2pi;
     }
     if(is_vbar && is_n) {
       if(npi==1) return kXSecTwkDial_RvbarnNC1pi;
       if(npi==2) return kXSecTwkDial_RvbarnNC2pi;
     }
   }//nc

   return kNullSystematic;
 }
 //......................................................................................

private:

  /// Map that defines conversions between GSyst_t values and their string
  /// representations
  static std::map<GSyst_t, std::string> fGSystToStringMap;

  /// Helper function to initialize the GSyst_t <-> std::string conversion map
  static std::map<GSyst_t, std::string> BuildGSystToStringMap();

};

} // rew   namespace
} // genie namespace

#endif
