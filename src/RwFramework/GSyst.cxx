//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// GENIE/Reweight includes
#include "RwFramework/GSyst.h"

using namespace genie;
using namespace genie::rew;

std::map<GSyst_t, std::string> GSyst::fGSystToStringMap
  = GSyst::BuildGSystToStringMap();

// There are better ways of doing this in C++11. I'd recommend
// revisiting this for GENIE 4. - S. Gardiner
std::map<GSyst_t, std::string> GSyst::BuildGSystToStringMap() {
  std::map<GSyst_t, std::string> temp_map;

  temp_map[ kXSecTwkDial_MaNCEL ]           = "MaNCEL";
  temp_map[ kXSecTwkDial_EtaNCEL ]          = "EtaNCEL";
  temp_map[ kXSecTwkDial_NormCCQE ]         = "NormCCQE";
  temp_map[ kXSecTwkDial_NormCCQEenu ]      = "NormCCQEenu";
  temp_map[ kXSecTwkDial_MaCCQE ]           = "MaCCQE";
  temp_map[ kXSecTwkDial_MaCCQEshape ]      = "MaCCQEshape";
  temp_map[ kXSecTwkDial_E0CCQE ]           = "E0CCQE";
  temp_map[ kXSecTwkDial_E0CCQEshape ]      = "E0CCQEshape";
  temp_map[ kXSecTwkDial_ZNormCCQE ]        = "ZNormCCQE";
  temp_map[ kXSecTwkDial_ZExpA1CCQE ]       = "ZExpA1CCQE";
  temp_map[ kXSecTwkDial_ZExpA2CCQE ]       = "ZExpA2CCQE";
  temp_map[ kXSecTwkDial_ZExpA3CCQE ]       = "ZExpA3CCQE";
  temp_map[ kXSecTwkDial_ZExpA4CCQE ]       = "ZExpA4CCQE";
  temp_map[ kXSecTwkDial_AxFFCCQEshape ]    = "AxFFCCQEshape";
  temp_map[ kXSecTwkDial_VecFFCCQEshape ]   = "VecFFCCQEshape";
  temp_map[ kXSecTwkDial_NormCCRES ]        = "NormCCRES";
  temp_map[ kXSecTwkDial_MaCCRESshape ]     = "MaCCRESshape";
  temp_map[ kXSecTwkDial_MvCCRESshape ]     = "MvCCRESshape";
  temp_map[ kXSecTwkDial_MaCCRES ]          = "MaCCRES";
  temp_map[ kXSecTwkDial_MvCCRES ]          = "MvCCRES";
  temp_map[ kXSecTwkDial_NormNCRES ]        = "NormNCRES";
  temp_map[ kXSecTwkDial_MaNCRESshape ]     = "MaNCRESshape";
  temp_map[ kXSecTwkDial_MvNCRESshape ]     = "MvNCRESshape";
  temp_map[ kXSecTwkDial_MaNCRES ]          = "MaNCRES";
  temp_map[ kXSecTwkDial_MvNCRES ]          = "MvNCRES";
  temp_map[ kXSecTwkDial_MaCOHpi ]          = "MaCOHpi";
  temp_map[ kXSecTwkDial_R0COHpi ]          = "R0COHpi";
  temp_map[ kXSecTwkDial_RvpCC1pi ]         = "NonRESBGvpCC1pi";
  temp_map[ kXSecTwkDial_RvpCC2pi ]         = "NonRESBGvpCC2pi";
  temp_map[ kXSecTwkDial_RvpNC1pi ]         = "NonRESBGvpNC1pi";
  temp_map[ kXSecTwkDial_RvpNC2pi ]         = "NonRESBGvpNC2pi";
  temp_map[ kXSecTwkDial_RvnCC1pi ]         = "NonRESBGvnCC1pi";
  temp_map[ kXSecTwkDial_RvnCC2pi ]         = "NonRESBGvnCC2pi";
  temp_map[ kXSecTwkDial_RvnNC1pi ]         = "NonRESBGvnNC1pi";
  temp_map[ kXSecTwkDial_RvnNC2pi ]         = "NonRESBGvnNC2pi";
  temp_map[ kXSecTwkDial_RvbarpCC1pi ]      = "NonRESBGvbarpCC1pi";
  temp_map[ kXSecTwkDial_RvbarpCC2pi ]      = "NonRESBGvbarpCC2pi";
  temp_map[ kXSecTwkDial_RvbarpNC1pi ]      = "NonRESBGvbarpNC1pi";
  temp_map[ kXSecTwkDial_RvbarpNC2pi ]      = "NonRESBGvbarpNC2pi";
  temp_map[ kXSecTwkDial_RvbarnCC1pi ]      = "NonRESBGvbarnCC1pi";
  temp_map[ kXSecTwkDial_RvbarnCC2pi ]      = "NonRESBGvbarnCC2pi";
  temp_map[ kXSecTwkDial_RvbarnNC1pi ]      = "NonRESBGvbarnNC1pi";
  temp_map[ kXSecTwkDial_RvbarnNC2pi ]      = "NonRESBGvbarnNC2pi";
  temp_map[ kXSecTwkDial_AhtBY ]            = "AhtBY";
  temp_map[ kXSecTwkDial_BhtBY ]            = "BhtBY";
  temp_map[ kXSecTwkDial_CV1uBY ]           = "CV1uBY";
  temp_map[ kXSecTwkDial_CV2uBY ]           = "CV2uBY";
  temp_map[ kXSecTwkDial_AhtBYshape ]       = "AhtBYshape";
  temp_map[ kXSecTwkDial_BhtBYshape ]       = "BhtBYshape";
  temp_map[ kXSecTwkDial_CV1uBYshape ]      = "CV1uBYshape";
  temp_map[ kXSecTwkDial_CV2uBYshape ]      = "CV2uBYshape";
  temp_map[ kXSecTwkDial_NormDISCC ]        = "NormDISCC";
  temp_map[ kXSecTwkDial_RnubarnuCC ]       = "RnubarnuCC";
  temp_map[ kXSecTwkDial_DISNuclMod ]       = "DISNuclMod";
  temp_map[ kXSecTwkDial_NC ]               = "NC";
  temp_map[ kHadrAGKYTwkDial_xF1pi ]        = "AGKYxF1pi";
  temp_map[ kHadrAGKYTwkDial_pT1pi ]        = "AGKYpT1pi";
  temp_map[ kHadrNuclTwkDial_FormZone ]     = "FormZone";
  temp_map[ kINukeTwkDial_MFP_pi ]          = "MFP_pi";
  temp_map[ kINukeTwkDial_MFP_N ]           = "MFP_N";
  temp_map[ kINukeTwkDial_FrCEx_pi ]        = "FrCEx_pi";
//temp_map[ kINukeTwkDial_FrElas_pi ]       = "FrElas_pi";
  temp_map[ kINukeTwkDial_FrInel_pi ]       = "FrInel_pi";
  temp_map[ kINukeTwkDial_FrAbs_pi ]        = "FrAbs_pi";
  temp_map[ kINukeTwkDial_FrPiProd_pi ]     = "FrPiProd_pi";
  temp_map[ kINukeTwkDial_FrCEx_N ]         = "FrCEx_N";
//temp_map[ kINukeTwkDial_FrElas_N ]        = "FrElas_N";
  temp_map[ kINukeTwkDial_FrInel_N ]        = "FrInel_N";
  temp_map[ kINukeTwkDial_FrAbs_N ]         = "FrAbs_N";
  temp_map[ kINukeTwkDial_FrPiProd_N ]      = "FrPiProd_N";
  temp_map[ kSystNucl_CCQEPauliSupViaKF ]   = "CCQEPauliSupViaKF";
  temp_map[ kSystNucl_CCQEMomDistroFGtoSF ] = "CCQEMomDistroFGtoSF";
  temp_map[ kRDcyTwkDial_BR1gamma ]         = "RDecBR1gamma";
  temp_map[ kRDcyTwkDial_BR1eta ]           = "RDecBR1eta";
  temp_map[ kRDcyTwkDial_Theta_Delta2Npi ]  = "Theta_Delta2Npi";
  temp_map[ kXSecTwkDial_EmpMEC_Mq2d ]      = "EmpMEC_Mq2d";
  temp_map[ kXSecTwkDial_EmpMEC_Mass ]      = "EmpMEC_Mass";
  temp_map[ kXSecTwkDial_EmpMEC_Width ]     = "EmpMEC_Width";
  temp_map[ kXSecTwkDial_EmpMEC_FracPN_NC ] = "EmpMEC_FracPN_NC";
  temp_map[ kXSecTwkDial_EmpMEC_FracPN_CC ] = "EmpMEC_FracPN_CC";
  temp_map[ kXSecTwkDial_EmpMEC_FracCCQE ]  = "EmpMEC_FracCCQE";
  temp_map[ kXSecTwkDial_EmpMEC_FracNCQE ]  = "EmpMEC_FracNCQE";
  temp_map[ kXSecTwkDial_EmpMEC_FracPN_EM ] = "EmpMEC_FracPN_EM";
  temp_map[ kXSecTwkDial_EmpMEC_FracEMQE ]  = "EmpMEC_FracEMQE";

  return temp_map;
}
