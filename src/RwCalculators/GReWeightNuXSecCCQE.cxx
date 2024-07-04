//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TFile.h>
#include <TNtupleD.h>
#include <cstdlib>
#include <sstream>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Registry/Registry.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightNuXSecCCQE.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;
using std::ostringstream;

static const char* kModelDipole        = "genie::DipoleAxialFormFactorModel";
static const char* kModelZExp          = "genie::ZExpAxialFormFactorModel";
static const char* kModelMArunAxial    = "genie::MArunAxialFormFactorModel";

const int GReWeightNuXSecCCQE::kModeMa;
const int GReWeightNuXSecCCQE::kModeNormAndMaShape;
const int GReWeightNuXSecCCQE::kModeZExp;
const int GReWeightNuXSecCCQE::fZExpMaxSyst;

//_______________________________________________________________________________________
GReWeightNuXSecCCQE::GReWeightNuXSecCCQE() :
GReWeightModel("CCQE"),
fManualModelName(),
fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQE::GReWeightNuXSecCCQE(std::string model, std::string type) :
GReWeightModel("CCQE"),
fManualModelName(model),
fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQE::~GReWeightNuXSecCCQE()
{
  if ( fXSecModelConfig ) delete fXSecModelConfig;

  if ( fXSecModel ) delete fXSecModel;
  if ( fXSecModelDef ) delete fXSecModelDef;
  if ( fXSecModelDefNoRPA ) delete fXSecModelDefNoRPA;
  if ( fXSecModelDefCoulomb ) delete fXSecModelDefCoulomb;

  if ( fXSecIntegrator ) delete fXSecIntegrator;
  if ( fXSecIntegratorDef ) delete fXSecIntegratorDef;

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQE::IsHandled(GSyst_t syst) const
{
// read form factor model and compare to mode
   bool handle;

   switch(syst) {

     case ( kXSecTwkDial_RPA_CCQE ) :
     case ( kXSecTwkDial_CoulombCCQE ) :
       handle = true;
       break;

     case ( kXSecTwkDial_NormCCQE    ) :
       if(fMode==kModeNormAndMaShape && (fModelIsDipole || fModelIsMArunAxial))
       {
          handle = true;
       } else {
          handle = false;
       }
       break;

     case ( kXSecTwkDial_MaCCQEshape ) :
     case ( kXSecTwkDial_E0CCQEshape ) :
       if(fMode==kModeNormAndMaShape && (fModelIsDipole || fModelIsMArunAxial) )
       {
          handle = true;
       } else {
          handle = false;
       }
       break;

     case ( kXSecTwkDial_MaCCQE ) :
     case ( kXSecTwkDial_E0CCQE ) :
       if(fMode==kModeMa && (fModelIsDipole || fModelIsMArunAxial))
       {
          handle = true;
       } else {
          handle = false;
       }
       break;

     case ( kXSecTwkDial_ZNormCCQE    ) :
       if(fMode==kModeZExp && fModelIsZExp)
       {
          handle = true;
       } else {
          handle = false;
       }
       break;

     case ( kXSecTwkDial_ZExpA1CCQE ):
     case ( kXSecTwkDial_ZExpA2CCQE ):
     case ( kXSecTwkDial_ZExpA3CCQE ):
     case ( kXSecTwkDial_ZExpA4CCQE ):
       if(fMode==kModeZExp && fModelIsZExp)
       {
          handle = true;
       } else {
          handle = false;
       }
       break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQE::AppliesTo(const EventRecord &event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if (type==kScQuasiElastic && is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst))
  {
    LOG("ReW",pWARN) << "Systematic " << GSyst::AsString(syst) << " is not handled for algorithm "
      << fFFModel << " and mode " << fMode;
    return;
  }

    switch(syst) {
      case ( kXSecTwkDial_NormCCQE ) :
      case ( kXSecTwkDial_ZNormCCQE ) :
        fNormTwkDial = twk_dial;
        break;
      case ( kXSecTwkDial_MaCCQEshape ) :
      case ( kXSecTwkDial_MaCCQE ) :
        fMaTwkDial = twk_dial;
        break;
      case ( kXSecTwkDial_E0CCQEshape ) :
      case ( kXSecTwkDial_E0CCQE ) :
        fE0TwkDial = twk_dial;
        break;
      case ( kXSecTwkDial_ZExpA1CCQE ) :
        if(fZExpMaxCoef>0){ fZExpTwkDial[0] = twk_dial; }
        break;
      case ( kXSecTwkDial_ZExpA2CCQE ) :
        if(fZExpMaxCoef>1){ fZExpTwkDial[1] = twk_dial; }
        break;
      case ( kXSecTwkDial_ZExpA3CCQE ) :
        if(fZExpMaxCoef>2){ fZExpTwkDial[2] = twk_dial; }
        break;
      case ( kXSecTwkDial_ZExpA4CCQE ) :
        if(fZExpMaxCoef>3){ fZExpTwkDial[3] = twk_dial; }
        break;
      case ( kXSecTwkDial_RPA_CCQE ) :
        fRPATwkDial = twk_dial;
        break;
      case ( kXSecTwkDial_CoulombCCQE ) :
        fCoulombTwkDial = twk_dial;
        break;
      default:
        break;
    }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.;
  fMaCurr      = fMaDef;
  fE0Curr      = fE0Def;

  for (int i=0;i<fZExpMaxSyst;i++)
  {
    fZExpTwkDial[i] = 0.;
    fZExpCurr   [i] = fZExpDef[i];
  }

  fRPATwkDial = 0.;
  fCoulombTwkDial = 0.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  if(fMode==kModeMa && fModelIsDipole) {
     int    sign_matwk = utils::rew::Sign(fMaTwkDial);
     double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQE, sign_matwk);
     fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
     fMaCurr = TMath::Max(0., fMaCurr  );
  }
  else
  if(fMode==kModeNormAndMaShape && fModelIsDipole) {
     int    sign_normtwk = utils::rew::Sign(fNormTwkDial);
     int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormCCQE,  sign_normtwk);
     double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQEshape, sign_mashtwk);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
     fMaCurr   = fMaDef   * (1. + fMaTwkDial   * fracerr_mash);
     fNormCurr = TMath::Max(0., fNormCurr);
     fMaCurr   = TMath::Max(0., fMaCurr  );
  }
  else
  if(fMode==kModeMa && fModelIsMArunAxial) {
     int    sign_matwk = utils::rew::Sign(fMaTwkDial);
     int    sign_e0twk = utils::rew::Sign(fE0TwkDial);
     double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQE, sign_matwk);
     double fracerr_e0 = fracerr->OneSigmaErr(kXSecTwkDial_E0CCQE, sign_e0twk);
     fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
     fMaCurr = TMath::Max(0., fMaCurr  );
     fE0Curr = fE0Def * (1. + fE0TwkDial * fracerr_e0);
     fE0Curr = TMath::Max(0., fE0Curr  );
  }
  else
  if(fMode==kModeNormAndMaShape && fModelIsMArunAxial) {
     int    sign_normtwk = utils::rew::Sign(fNormTwkDial);
     int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
     int    sign_e0shtwk = utils::rew::Sign(fE0TwkDial);
     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormCCQE,  sign_normtwk);
     double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQEshape, sign_mashtwk);
     double fracerr_e0sh = fracerr->OneSigmaErr(kXSecTwkDial_E0CCQEshape, sign_e0shtwk);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
     fMaCurr   = fMaDef   * (1. + fMaTwkDial   * fracerr_mash);
     fE0Curr = fE0Def * (1. + fE0TwkDial * fracerr_e0sh);
     fNormCurr = TMath::Max(0., fNormCurr);
     fMaCurr   = TMath::Max(0., fMaCurr  );
     fE0Curr = TMath::Max(0., fE0Curr  );
  }
  else
  if(fMode==kModeZExp && fModelIsZExp) {
     int     sign_twk = 0;
     int     sign_normtwk = utils::rew::Sign(fNormTwkDial);
     double  fracerr_zexp = 0.;
     double  fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_ZNormCCQE, sign_normtwk);
     GSyst_t syst;
     // loop over all indices and update each
     for (int i=0;i<fZExpMaxCoef;i++)
     {
       switch(i){
         case 0: syst = kXSecTwkDial_ZExpA1CCQE; break;
         case 1: syst = kXSecTwkDial_ZExpA2CCQE; break;
         case 2: syst = kXSecTwkDial_ZExpA3CCQE; break;
         case 3: syst = kXSecTwkDial_ZExpA4CCQE; break;
         default: return; break;
       }
       sign_twk = utils::rew::Sign(fZExpTwkDial[i]);
       fracerr_zexp = fracerr->OneSigmaErr(syst, sign_twk);
       fZExpCurr[i] = fZExpDef[i] * (1. + fZExpTwkDial[i] * fracerr_zexp);
     }
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
     fNormCurr = TMath::Max(0., fNormCurr);
  }
  else {
    return;
  }

  Registry r("GReWeightNuXSecCCQE",false);
  //~ Registry r(fXSecModel->GetConfig());
  if (fMode==kModeMa || fMode==kModeNormAndMaShape)
  {
    r.Set(fMaPath, fMaCurr);
    r.Set(fE0Path, fE0Curr);
  }
  else
  if (fMode==kModeZExp)
  {
    ostringstream alg_key;
    for (int i=0;i<fZExpMaxCoef;i++)
    {
      alg_key.str(""); // algorithm key for each coefficient
      alg_key << fZExpPath << "QEL-Z_A" << i+1;
      r.Set(alg_key.str(), fZExpCurr[i]);
    }
  }
  fXSecModel->Configure(r);

  // Update the tweaked value of the Coulomb scaling factor
  int sign_twk_C = utils::rew::Sign( fCoulombTwkDial );
  double frac_sigma_C = fracerr->OneSigmaErr( kXSecTwkDial_CoulombCCQE,
    sign_twk_C );
  double scaleC = (1. + fCoulombTwkDial * frac_sigma_C);

  // Update the configuration for the Coulomb-tweaked cross section model
  Registry coulomb_reg = fXSecModelDef->GetConfig();
  coulomb_reg.Set( "CoulombScale", scaleC );
  fXSecModelDefCoulomb->Configure( coulomb_reg );
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeight(const genie::EventRecord & event)
{
  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if ( !is_qe || !is_cc ) return 1.;

  bool charm = event.Summary()->ExclTag().IsCharmEvent(); // skip CCQE charm
  if ( charm ) return 1.;

  // Skip CCQE strange
  bool strange = event.Summary()->ExclTag().IsStrangeEvent();
  if ( strange ) return 1.;

  // Skip any other CCQE channels that do not produce a final-state nucleon
  int final_nucleon_pdgc = event.Summary()->RecoilNucleonPdg();
  if ( !pdg::IsNucleon(final_nucleon_pdgc) ) return 1.;

  int nupdg = event.Probe()->Pdg();
  if ( nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if ( nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if ( nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if ( nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  double wght = this->CalcWeightRPA( event );

  wght *= this->CalcWeightCoulomb( event );

  if ( fMode==kModeMa && (fModelIsDipole || fModelIsMArunAxial) ) {
     wght *= this->CalcWeightMa(event);
     return wght;
  }
  else
  if ( fMode==kModeNormAndMaShape && (fModelIsDipole || fModelIsMArunAxial) ) {
     wght *= this->CalcWeightNorm( event ) * this->CalcWeightMaShape( event );
     return wght;
  }
  else
  if ( fMode==kModeZExp && fModelIsZExp ) {
     wght *= this->CalcWeightNorm( event ) * this->CalcWeightZExp( event );
     return wght;
  }

  return 1.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/QEL-CC");

  AlgId id(xsec_alg);

  AlgId twk_id(id);
  if (fManualModelName.size()) {
    twk_id = AlgId(fManualModelName,fManualModelType);
  }

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * alg_def = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(alg_def);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg_twk = algf->AdoptAlgorithm(twk_id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg_twk);
  fXSecModel->AdoptSubstructure();

  // Get the Algorithm objects that should be used to integrate the cross
  // sections. Use the "ReweightShape" configuration, which turns off averaging
  // over the initial state nucleon distribution.
  AlgId alg_def_integ_ID( alg_def->GetConfig().GetAlg("XSec-Integrator").name,
    "ReweightShape");

  fXSecIntegratorDef = dynamic_cast<XSecIntegratorI*>(
    algf->AdoptAlgorithm(alg_def_integ_ID));

  assert( fXSecIntegratorDef );

  AlgId alg_twk_integ_ID( alg_twk->GetConfig().GetAlg("XSec-Integrator").name,
    "ReweightShape");

  fXSecIntegrator = dynamic_cast<XSecIntegratorI*>(
    algf->AdoptAlgorithm(alg_twk_integ_ID));

  assert( fXSecIntegrator );

  // Check what kind of form factors we're using in the tweaked cross section
  // model
  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
  fFFModel = fXSecModelConfig->GetAlg("FormFactorsAlg/AxialFormFactorModel").name;

  fModelIsDipole    = (strcmp(fFFModel.c_str(),kModelDipole) == 0);
  fModelIsZExp      = (strcmp(fFFModel.c_str(),kModelZExp  ) == 0);
  fModelIsMArunAxial = (strcmp(fFFModel.c_str(),kModelMArunAxial  ) == 0);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  //this->SetMaPath("FormFactorsAlg/Ma");
  this->SetMaPath("FormFactorsAlg/AxialFormFactorModel/QEL-Ma");
  this->SetE0Path("FormFactorsAlg/AxialFormFactorModel/QEL-E0");
  this->SetZExpPath("FormFactorsAlg/AxialFormFactorModel/");

  if (fModelIsDipole)
  {
    this->SetMode(kModeNormAndMaShape);
    fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
    fZExpMaxCoef = 0;
  }
  else if (fModelIsMArunAxial)
  {
    this->SetMode(kModeNormAndMaShape);
    fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
    fE0Def       = fXSecModelConfig->GetDouble(fE0Path);
    fZExpMaxCoef = 0;
  }
  else if (fModelIsZExp)
  {
    this->SetMode(kModeZExp);
    fMaDef = 0.;
    fZExpMaxCoef = TMath::Min(fXSecModelConfig->GetInt(fZExpPath + "QEL-Kmax"),
    this->fZExpMaxSyst);
  }

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.;
  fMaCurr      = fMaDef;
  fE0TwkDial   = 0.;
  fE0Curr      = fE0Def;

  ostringstream alg_key;
  for (int i=0;i<fZExpMaxSyst;i++)
  {
    alg_key.str("");
    alg_key << fZExpPath << "QEL-Z_A" << i+1;
    if (fModelIsZExp && i < fZExpMaxCoef)
    { fZExpDef[i] = fXSecModelConfig->GetDouble(alg_key.str()); }
    else
    { fZExpDef[i] = 0.; }
    fZExpTwkDial[i] = 0.;
    fZExpCurr   [i] = fZExpDef[i];
  }

  fRPATwkDial = 0.;
  fCoulombTwkDial = 0.;

  // Make a copy of the default cross section model that has RPA turned off
  Algorithm* alg_def_copy = algf->AdoptAlgorithm( id );
  fXSecModelDefNoRPA = dynamic_cast<XSecAlgorithmI*>( alg_def_copy );
  fXSecModelDefNoRPA->AdoptSubstructure();

  Registry rpa_off_reg = fXSecModelDef->GetConfig();
  rpa_off_reg.Set("RPA", false);

  fXSecModelDefNoRPA->Configure( rpa_off_reg );

  // Make a copy of the default cross section model with a tweaked
  // value of the scaling factor for the Coulomb potential
  Algorithm* alg_def_copy2 = algf->AdoptAlgorithm( id );
  fXSecModelDefCoulomb = dynamic_cast<XSecAlgorithmI*>( alg_def_copy2 );
  fXSecModelDefCoulomb->AdoptSubstructure();

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  fTestFile = new TFile("./ccqe_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght");
#endif
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightNorm(const genie::EventRecord & /*event*/)
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightMa(const genie::EventRecord & event)
{
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum) || ((TMath::Abs(fE0TwkDial) > controls::kASmallNum) && fModelIsMArunAxial);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  // Retrieve the kinematic phase space used to generate the event
  const KinePhaseSpace_t phase_space = event.DiffXSecVars();
  double old_xsec = event.DiffXSec();

  if (phase_space == kPSQ2fE) {
    interaction->SetBit(kIAssumeFreeNucleon);
  }

  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
    if (fNWeightChecksDone < fNWeightChecksToDo) {
      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
        fFailedWeightCheck = true;
      }
      fNWeightChecksDone++;
    }
    if(!fUseOldWeightFromFile) {
      old_xsec = calc_old_xsec;
    }
  }

  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightMaShape(const genie::EventRecord & event)
{
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum) || ((TMath::Abs(fE0TwkDial) > controls::kASmallNum) && fModelIsMArunAxial);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  // Retrieve the kinematic phase space used to generate the event
  const KinePhaseSpace_t phase_space = event.DiffXSecVars();
  double old_xsec = event.DiffXSec();

  if (phase_space == kPSQ2fE) {
    interaction->SetBit(kIAssumeFreeNucleon);
  }

  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
    if (fNWeightChecksDone < fNWeightChecksToDo) {
      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
      }
      fNWeightChecksDone++;
    }
    if(!fUseOldWeightFromFile) {
      old_xsec = calc_old_xsec;
    }
  }
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

  double old_integrated_xsec = fXSecIntegratorDef->Integrate(fXSecModelDef, interaction);
  double new_integrated_xsec = fXSecIntegrator->Integrate(fXSecModel, interaction);
  assert(new_integrated_xsec > 0);
  new_weight *= (old_integrated_xsec/new_integrated_xsec);

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " << old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (new) = " << new_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(E,Q2,new_weight);
#endif

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightZExp(const genie::EventRecord & event)
{
  // very similar to CalcWeightMa
  bool tweaked = false;
  for (int i=0;i<fZExpMaxCoef;i++)
  {
    tweaked = tweaked || (TMath::Abs(fZExpTwkDial[i]) > controls::kASmallNum);
  }
  if(!tweaked) { return 1.0; }

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  // Retrieve the kinematic phase space used to generate the event
  const KinePhaseSpace_t phase_space = event.DiffXSecVars();
  double old_xsec = event.DiffXSec();

  if (phase_space == kPSQ2fE) {
    interaction->SetBit(kIAssumeFreeNucleon);
  }

  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
    if (fNWeightChecksDone < fNWeightChecksToDo) {
      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
      }
      fNWeightChecksDone++;
    }
    if(!fUseOldWeightFromFile) {
      old_xsec = calc_old_xsec;
    }
  }
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  double new_weight = old_weight * (new_xsec/old_xsec);

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightRPA(const genie::EventRecord& event)
{
  // Only bother to do the calculation if the tweak dial value was noticeably
  // altered
  bool tweaked = ( std::abs(fRPATwkDial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  Interaction* interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  // Retrieve the kinematic phase space used to generate the event
  const KinePhaseSpace_t phase_space = event.DiffXSecVars();
  double old_xsec = event.DiffXSec();

  if ( phase_space == kPSQ2fE ) interaction->SetBit(kIAssumeFreeNucleon);

  if  ( !fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo ) {
    double calc_old_xsec = fXSecModelDef->XSec( interaction, phase_space );
    if ( fNWeightChecksDone < fNWeightChecksToDo ) {
      if ( std::abs(calc_old_xsec - old_xsec) / old_xsec > controls::kASmallNum ) {
        LOG("ReW", pWARN) << "Warning - default dxsec does not match dxsec"
          << " saved in tree. Does the config match?";
      }
      ++fNWeightChecksDone;
    }
    if ( !fUseOldWeightFromFile ) old_xsec = calc_old_xsec;
  }

  assert( old_xsec > 0. );

  double rpa_off_xsec = fXSecModelDefNoRPA->XSec( interaction, phase_space );

  LOG("ReW", pDEBUG) << "old_xsec = " << old_xsec << ", RPA off = " << rpa_off_xsec;

  double new_xsec = fRPATwkDial*rpa_off_xsec + (1. - fRPATwkDial)*old_xsec;
  double weight = new_xsec / old_xsec;

  // Allow extrapolations outside of the nominal [0, 1] range of the RPA_CCQE
  // dial by zeroing out negative weights
  if ( weight < 0. ) weight = 0.;

  interaction->KinePtr()->ClearRunningValues();

  if ( phase_space == kPSQ2fE ) interaction->ResetBit(kIAssumeFreeNucleon);

  return weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightCoulomb(const genie::EventRecord& event)
{
  // Only bother to do the calculation if the tweak dial value was noticeably
  // altered
  bool tweaked = ( std::abs(fCoulombTwkDial) > controls::kASmallNum );
  if ( !tweaked ) return 1.;

  Interaction* interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  // Retrieve the kinematic phase space used to generate the event
  const KinePhaseSpace_t phase_space = event.DiffXSecVars();
  double old_xsec = event.DiffXSec();

  if ( phase_space == kPSQ2fE ) interaction->SetBit(kIAssumeFreeNucleon);

  if  ( !fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo ) {
    double calc_old_xsec = fXSecModelDef->XSec( interaction, phase_space );
    if ( fNWeightChecksDone < fNWeightChecksToDo ) {
      if ( std::abs(calc_old_xsec - old_xsec) / old_xsec > controls::kASmallNum ) {
        LOG("ReW", pWARN) << "Warning - default dxsec does not match dxsec"
          << " saved in tree. Does the config match?";
      }
      ++fNWeightChecksDone;
    }
    if ( !fUseOldWeightFromFile ) old_xsec = calc_old_xsec;
  }

  assert( old_xsec > 0. );

  double new_xsec = fXSecModelDefCoulomb->XSec( interaction, phase_space );

  LOG("ReW", pDEBUG) << "old_xsec = " << old_xsec << ", tweaked Coulomb"
    " xsec = " << new_xsec;

  double weight = new_xsec / old_xsec;

  interaction->KinePtr()->ClearRunningValues();

  if ( phase_space == kPSQ2fE ) interaction->ResetBit(kIAssumeFreeNucleon);

  return weight;
}
