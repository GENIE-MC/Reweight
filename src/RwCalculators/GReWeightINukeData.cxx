//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>

#include <TTree.h>
#include <TSystem.h>

#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "RwFramework/GSyst.h"
#include "RwCalculators/GReWeightINukeData.h"

using namespace genie;
using namespace genie::rew;

// the Instance
GReWeightINukeData *theInstance = 0;

const GReWeightINukeData *GReWeightINukeData::Instance() {

  if (theInstance == 0) theInstance = new GReWeightINukeData;

  return theInstance;

}

GReWeightINukeData::GReWeightINukeData() {
  // Load location of hadron data
  std::string data_dir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
             string(gSystem->Getenv("GINUKEHADRONDATA")) :
             string(gSystem->Getenv("GENIE")) + string("/data/evgen/intranuke");


  // data files with variations
  std::string datafile_NAG4 = data_dir + "/tot_xsec/intranuke-fractions-NA2025G4.dat";
  std::string datafile_NAINCL = data_dir + "/tot_xsec/intranuke-fractions-NA2025INCL.dat";

  // verify files exist
  assert( ! gSystem->AccessPathName(datafile_NAG4.  c_str()) );
  assert( ! gSystem->AccessPathName(datafile_NAINCL.  c_str()) );

  // load into TTrees
  TTree data_NAG4;
  TTree data_NAINCL;

  data_NAG4.ReadFile(datafile_NAG4.c_str(), 
                   "ke/D:pA_tot/D:pA_el/D:pA_inel/D:pA_cex/D:pA_abs/D:pA_pipro/D:pA_xsec/D");
  data_NAINCL.ReadFile(datafile_NAINCL.c_str(), 
                   "ke/D:pA_tot/D:pA_el/D:pA_inel/D:pA_cex/D:pA_abs/D:pA_pipro/D:pA_xsec/D");

  // Load splines
  fG4PA_Tot = new Spline(&data_NAG4, "ke:pA_xsec");
  fG4FracPA_Tot = new Spline(&data_NAG4, "ke:pA_tot");
  fG4FracPA_Inel = new Spline(&data_NAG4, "ke:pA_inel");
  fG4FracPA_CEx = new Spline(&data_NAG4, "ke:pA_cex");
  fG4FracPA_Abs = new Spline(&data_NAG4, "ke:pA_abs");
  fG4FracPA_PiPro = new Spline(&data_NAG4, "ke:pA_pipro");

  fG4NA_Tot = new Spline(&data_NAG4, "ke:pA_xsec");
  fG4FracNA_Tot = new Spline(&data_NAG4, "ke:pA_tot");
  fG4FracNA_Inel = new Spline(&data_NAG4, "ke:pA_inel");
  fG4FracNA_CEx = new Spline(&data_NAG4, "ke:pA_cex");
  fG4FracNA_Abs = new Spline(&data_NAG4, "ke:pA_abs");
  fG4FracNA_PiPro = new Spline(&data_NAG4, "ke:pA_pipro");

  fINCLPA_Tot = new Spline(&data_NAINCL, "ke:pA_xsec");
  fINCLFracPA_Tot = new Spline(&data_NAINCL, "ke:pA_tot");
  fINCLFracPA_Inel = new Spline(&data_NAINCL, "ke:pA_inel");
  fINCLFracPA_CEx = new Spline(&data_NAINCL, "ke:pA_cex");
  fINCLFracPA_Abs = new Spline(&data_NAINCL, "ke:pA_abs");
  fINCLFracPA_PiPro = new Spline(&data_NAINCL, "ke:pA_pipro");

  fINCLNA_Tot = new Spline(&data_NAINCL, "ke:pA_xsec");
  fINCLFracNA_Tot = new Spline(&data_NAINCL, "ke:pA_tot");
  fINCLFracNA_Inel = new Spline(&data_NAINCL, "ke:pA_inel");
  fINCLFracNA_CEx = new Spline(&data_NAINCL, "ke:pA_cex");
  fINCLFracNA_Abs = new Spline(&data_NAINCL, "ke:pA_abs");
  fINCLFracNA_PiPro = new Spline(&data_NAINCL, "ke:pA_pipro");

}

// cleanup splines
GReWeightINukeData::~GReWeightINukeData() {
   delete fG4PA_Tot;
   delete fG4FracPA_Tot;
   delete fG4FracPA_Inel;
   delete fG4FracPA_CEx;
   delete fG4FracPA_Abs;
   delete fG4FracPA_PiPro;

   delete fG4NA_Tot;
   delete fG4FracNA_Tot;
   delete fG4FracNA_Inel;
   delete fG4FracNA_CEx;
   delete fG4FracNA_Abs;
   delete fG4FracNA_PiPro;

   delete fINCLPA_Tot;
   delete fINCLFracPA_Tot;
   delete fINCLFracPA_Inel;
   delete fINCLFracPA_CEx;
   delete fINCLFracPA_Abs;
   delete fINCLFracPA_PiPro;

   delete fINCLNA_Tot;
   delete fINCLFracNA_Tot;
   delete fINCLFracNA_Inel;
   delete fINCLFracNA_CEx;
   delete fINCLFracNA_Abs;
   delete fINCLFracNA_PiPro;
}
   
bool GReWeightINukeData::IsHandled(GSyst_t syst) const {
  switch(syst) {
    case ( kINukeTwkDial_FrCEx_N   ) :
    case ( kINukeTwkDial_FrInel_N  ) :
    case ( kINukeTwkDial_FrAbs_N   ) :
    case ( kINukeTwkDial_FrPiProd_N) :
      return true;
    default:
      return false;
  }
  return false;
}

//____________________________________________________________________________
double GReWeightINukeData::FateFraction(ModelSwitch_t model, GSyst_t syst, double kinE) const {
  double fate_frac = 0.0;

  // convert to MeV and
  double ke = kinE / units::MeV;
  ke = TMath::Max(GReWeightINukeData::fMinKinEnergy, ke);
  ke = TMath::Min(GReWeightINukeData::fMaxKinEnergy, ke);

  switch (syst) {
    case (kINukeTwkDial_FrCEx_N) :
    {
      if (model == kRwINukeG4) fate_frac = this->G4FracPA_CEx()->Evaluate(ke);
      else if (model == kRwINukeINCL) fate_frac = this->INCLFracPA_CEx()->Evaluate(ke);
    }
    break;

    case (kINukeTwkDial_FrInel_N) :
    {
      if (model == kRwINukeG4) fate_frac = this->G4FracPA_Inel()->Evaluate(ke);
      else if (model == kRwINukeINCL) fate_frac = this->INCLFracPA_Inel()->Evaluate(ke);
    }
    break;

    case (kINukeTwkDial_FrAbs_N) :
    {
      if (model == kRwINukeG4) fate_frac = this->G4FracPA_Abs()->Evaluate(ke);
      else if (model == kRwINukeINCL) fate_frac = this->INCLFracPA_Abs()->Evaluate(ke);
    }
    break;

    case (kINukeTwkDial_FrPiProd_N) :
    {
      if (model == kRwINukeG4) fate_frac = this->G4FracPA_PiPro()->Evaluate(ke);
      else if (model == kRwINukeINCL) fate_frac = this->INCLFracPA_PiPro()->Evaluate(ke);
    }
    break;

    default:
    {
      LOG("ReW", pDEBUG)
        << "Have reached default case and assigning fraction{fate} = 0";
      fate_frac = 0;
    }
    break;

  } // hadron_fate?

  return fate_frac;
}

