//____________________________________________________________________________
/*
   Copyright (c) 2003-2018, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org

Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
University of Liverpool & STFC Rutherford Appleton Lab

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
#include "GReWeightNuXSecCCQEELFF.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;
using std::ostringstream;

static const char* kModelZExp          = "genie::ZExpELFormFactorModel";

const int GReWeightNuXSecCCQEELFF::kModeZExp;

//_______________________________________________________________________________________
GReWeightNuXSecCCQEELFF::GReWeightNuXSecCCQEELFF() :
  GReWeightModel("CCQE"),
  fManualModelName(),
  fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEELFF::GReWeightNuXSecCCQEELFF(std::string model, std::string type) :
  GReWeightModel("CCQE"),
  fManualModelName(model),
  fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEELFF::~GReWeightNuXSecCCQEELFF()
{
  if ( fXSecModelConfig ) delete fXSecModelConfig;

  if ( fXSecModel ) delete fXSecModel;
  if ( fXSecModelDef ) delete fXSecModelDef;

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
bool GReWeightNuXSecCCQEELFF::IsHandled(GSyst_t syst) const
{
  // read form factor model and compare to mode
  bool handle;

  switch(syst) {
        // TODO add ZExp vector parameters
    case ( kXSecTwkDial_ZExpELFF_T0 ) :
    case ( kXSecTwkDial_ZExpELFF_Tcut ) :
    case ( kXSecTwkDial_ZExpELFF_Gep0 ) :
    case ( kXSecTwkDial_ZExpELFF_Gmp0 ) :
    case ( kXSecTwkDial_ZExpELFF_Gen0 ) :
    case ( kXSecTwkDial_ZExpELFF_Gmn0 ) :
    case ( kXSecTwkDial_ZExpELFF_AN1 ) :
    case ( kXSecTwkDial_ZExpELFF_AN2 ) :
    case ( kXSecTwkDial_ZExpELFF_AN3 ) :
    case ( kXSecTwkDial_ZExpELFF_AN4 ) :
    case ( kXSecTwkDial_ZExpELFF_AP1 ) :
    case ( kXSecTwkDial_ZExpELFF_AP2 ) :
    case ( kXSecTwkDial_ZExpELFF_AP3 ) :
    case ( kXSecTwkDial_ZExpELFF_AP4 ) :
    case ( kXSecTwkDial_ZExpELFF_BN1 ) :
    case ( kXSecTwkDial_ZExpELFF_BN2 ) :
    case ( kXSecTwkDial_ZExpELFF_BN3 ) :
    case ( kXSecTwkDial_ZExpELFF_BN4 ) :
    case ( kXSecTwkDial_ZExpELFF_BP1 ) :
    case ( kXSecTwkDial_ZExpELFF_BP2 ) :
    case ( kXSecTwkDial_ZExpELFF_BP3 ) :
    case ( kXSecTwkDial_ZExpELFF_BP4 ) : 
      if(fMode==kModeZExp && fModelIsZExp){
      handle = true;
      }else {
        handle = false;
      }
      break;
/*    case ( kXSecTwkDial_ZExpA1CCQE ):
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
*/
    default:
      handle = false;
      break;
  }

  return handle;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEELFF::AppliesTo(const EventRecord &event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if (type==kScQuasiElastic && is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEELFF::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst))
  {
    LOG("ReW",pWARN) << "Systematic " << GSyst::AsString(syst) << " is not handled for algorithm "
      << fFFModel << " and mode " << fMode;
    return;
  }
  switch(syst) {
    case ( kXSecTwkDial_ZExpELFF_T0 ) :
       fZExpParaTwkDial.fT0 = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_Tcut ) :
       fZExpParaTwkDial.fTcut = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_Gep0 ) :
       fZExpParaTwkDial.fGep0 = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_Gmp0 ) :
       fZExpParaTwkDial.fGmp0 = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_Gen0 ) :
       fZExpParaTwkDial.fGen0 = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_Gmn0 ) :
       fZExpParaTwkDial.fGmn0 = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AP1 ) :
       fZExpParaTwkDial.fZ_APn[0] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AP2 ) :
       fZExpParaTwkDial.fZ_APn[1] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AP3 ) :
       fZExpParaTwkDial.fZ_APn[2] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AP4 ) :
       fZExpParaTwkDial.fZ_APn[3] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AN1 ) :
       fZExpParaTwkDial.fZ_ANn[0] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AN2 ) :
       fZExpParaTwkDial.fZ_ANn[1] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AN3 ) :
       fZExpParaTwkDial.fZ_ANn[2] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_AN4 ) :
       fZExpParaTwkDial.fZ_ANn[3] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BP1 ) :
       fZExpParaTwkDial.fZ_BPn[0] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BP2 ) :
       fZExpParaTwkDial.fZ_BPn[1] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BP3 ) :
       fZExpParaTwkDial.fZ_BPn[2] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BP4 ) :
       fZExpParaTwkDial.fZ_BPn[3] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BN1 ) :
       fZExpParaTwkDial.fZ_BNn[0] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BN2 ) :
       fZExpParaTwkDial.fZ_BNn[1] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BN3 ) :
       fZExpParaTwkDial.fZ_BNn[2] = twk_dial; 
      break;
    case ( kXSecTwkDial_ZExpELFF_BN4 ) :
       fZExpParaTwkDial.fZ_BNn[3] = twk_dial; 
      break; 
//    case ( kXSecTwkDial_ZExpA2CCQE ) :
//      if(fZExpMaxCoef>1){ fZExpTwkDial[1] = twk_dial; }
//      break;
//    case ( kXSecTwkDial_ZExpA3CCQE ) :
//      if(fZExpMaxCoef>2){ fZExpTwkDial[2] = twk_dial; }
//      break;
//    case ( kXSecTwkDial_ZExpA4CCQE ) :
//      if(fZExpMaxCoef>3){ fZExpTwkDial[3] = twk_dial; }
//      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEELFF::Reset(void)
{
  fZExpPara.fQ4limit = fZExpParaDef.fQ4limit;
  fZExpPara.fKmax = fZExpParaDef.fKmax;
  fZExpPara.fT0 = fZExpParaDef.fT0;
  fZExpParaTwkDial.fT0 = 0.;
  fZExpPara.fTcut = fZExpParaDef.fTcut;
  fZExpParaTwkDial.fTcut = 0.;
  fZExpPara.fGep0 = fZExpParaDef.fGep0;
  fZExpParaTwkDial.fGep0 = 0.;
  fZExpPara.fGmp0 = fZExpParaDef.fGmp0;
  fZExpParaTwkDial.fGmp0 = 0.;
  fZExpPara.fGen0 = fZExpParaDef.fGen0;
  fZExpParaTwkDial.fGen0 = 0.;
  fZExpPara.fGmn0 = fZExpParaDef.fGmn0;
  fZExpParaTwkDial.fGmn0 = 0.;
  for(int i = 0; i < fZExpPara.fKmax; i++){
    fZExpPara.fZ_ANn[i] = fZExpParaDef.fZ_ANn[i];
    fZExpPara.fZ_APn[i] = fZExpParaDef.fZ_APn[i];
    fZExpPara.fZ_BNn[i] = fZExpParaDef.fZ_BNn[i];
    fZExpPara.fZ_BPn[i] = fZExpParaDef.fZ_BPn[i];
    fZExpParaTwkDial.fZ_ANn[i] = 0.;
    fZExpParaTwkDial.fZ_APn[i] = 0.;
    fZExpParaTwkDial.fZ_BNn[i] = 0.;
    fZExpParaTwkDial.fZ_APn[i] = 0.;
  }

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEELFF::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();
  if(fMode==kModeZExp && fModelIsZExp) {
    int     sign_twk = 0;
    double  fracerr_zexp = 0.;
//    double  fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_ZNormCCQE, sign_twk);
    sign_twk = utils::rew::Sign(fZExpParaTwkDial.fT0);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF_T0, sign_twk);
    fZExpPara.fT0 = fZExpParaDef.fT0 * (1. + fZExpParaTwkDial.fT0 * fracerr_zexp);
    sign_twk = utils::rew::Sign(fZExpParaTwkDial.fTcut);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF_Tcut, sign_twk);
    fZExpPara.fTcut = fZExpParaDef.fTcut * (1. + fZExpParaTwkDial.fTcut * fracerr_zexp);
    sign_twk = utils::rew::Sign(fZExpParaTwkDial.fGep0);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF_Gep0, sign_twk);
    fZExpPara.fGep0 = fZExpParaDef.fGep0 * (1. + fZExpParaTwkDial.fGep0 * fracerr_zexp);
    sign_twk = utils::rew::Sign(fZExpParaTwkDial.fGmp0);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF_Gmp0, sign_twk);
    fZExpPara.fGmp0 = fZExpParaDef.fGmp0 * (1. + fZExpParaTwkDial.fGmp0 * fracerr_zexp);
    sign_twk = utils::rew::Sign(fZExpParaTwkDial.fGen0);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF_Gen0, sign_twk);
    fZExpPara.fGen0 = fZExpParaDef.fGen0 * (1. + fZExpParaTwkDial.fGen0 * fracerr_zexp);
    sign_twk = utils::rew::Sign(fZExpParaTwkDial.fGmn0);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF_Gmn0, sign_twk);
    fZExpPara.fGmn0 = fZExpParaDef.fGmn0 * (1. + fZExpParaTwkDial.fGmn0 * fracerr_zexp);

    GSyst_t syst;
    // loop over all indices and update each
    for (int i=0;i<fZExpParaDef.fKmax;i++)
    {
      switch(i){
        case 0: syst = kXSecTwkDial_ZExpELFF_AP1; break;
        case 1: syst = kXSecTwkDial_ZExpELFF_AP2; break;
        case 2: syst = kXSecTwkDial_ZExpELFF_AP3; break;
        case 3: syst = kXSecTwkDial_ZExpELFF_AP4; break;
        default: return; break;
      }
      sign_twk = utils::rew::Sign(fZExpParaTwkDial.fZ_APn[i]);
      fracerr_zexp = fracerr->OneSigmaErr(syst, sign_twk);
      fZExpPara.fZ_APn[i] = fZExpParaDef.fZ_APn[i] * (1. + fZExpParaTwkDial.fZ_APn[i] * fracerr_zexp);
      switch(i){
        case 0: syst = kXSecTwkDial_ZExpELFF_BP1; break;
        case 1: syst = kXSecTwkDial_ZExpELFF_BP2; break;
        case 2: syst = kXSecTwkDial_ZExpELFF_BP3; break;
        case 3: syst = kXSecTwkDial_ZExpELFF_BP4; break;
        default: return; break;
      }
      sign_twk = utils::rew::Sign(fZExpParaTwkDial.fZ_BPn[i]);
      fracerr_zexp = fracerr->OneSigmaErr(syst, sign_twk);
      fZExpPara.fZ_BPn[i] = fZExpParaDef.fZ_BPn[i] * (1. + fZExpParaTwkDial.fZ_BPn[i] * fracerr_zexp);
      switch(i){
        case 0: syst = kXSecTwkDial_ZExpELFF_AN1; break;
        case 1: syst = kXSecTwkDial_ZExpELFF_AN2; break;
        case 2: syst = kXSecTwkDial_ZExpELFF_AN3; break;
        case 3: syst = kXSecTwkDial_ZExpELFF_AN4; break;
        default: return; break;
      }
      sign_twk = utils::rew::Sign(fZExpParaTwkDial.fZ_ANn[i]);
      fracerr_zexp = fracerr->OneSigmaErr(syst, sign_twk);
      fZExpPara.fZ_ANn[i] = fZExpParaDef.fZ_ANn[i] * (1. + fZExpParaTwkDial.fZ_ANn[i] * fracerr_zexp);
      switch(i){
        case 0: syst = kXSecTwkDial_ZExpELFF_BN1; break;
        case 1: syst = kXSecTwkDial_ZExpELFF_BN2; break;
        case 2: syst = kXSecTwkDial_ZExpELFF_BN3; break;
        case 3: syst = kXSecTwkDial_ZExpELFF_BN4; break;
        default: return; break;
      }
      sign_twk = utils::rew::Sign(fZExpParaTwkDial.fZ_BNn[i]);
      fracerr_zexp = fracerr->OneSigmaErr(syst, sign_twk);
      fZExpPara.fZ_BNn[i] = fZExpParaDef.fZ_BNn[i] * (1. + fZExpParaTwkDial.fZ_BNn[i] * fracerr_zexp);
    }
  }
  else {
    return;
  }


  Registry r("GReWeightNuXSecCCQEELFF",false);
  //~ Registry r(fXSecModel->GetConfig());
  if (fMode==kModeZExp)
  {
    ostringstream alg_key;
    alg_key.str(""); // algorithm key for each coefficient
    alg_key << fZExpPath << "QEL-T0";
    r.Set(alg_key.str(), fZExpPara.fT0);
    alg_key.str(""); // algorithm key for each coefficient
    alg_key << fZExpPath << "QEL-Tcut";
    r.Set(alg_key.str(), fZExpPara.fTcut);
    alg_key.str(""); // algorithm key for each coefficient
    alg_key << fZExpPath << "QEL-Gep0";
    r.Set(alg_key.str(), fZExpPara.fGep0);
    alg_key.str(""); // algorithm key for each coefficient
    alg_key << fZExpPath << "QEL-Gmp0";
    r.Set(alg_key.str(), fZExpPara.fGmp0);
    alg_key.str(""); // algorithm key for each coefficient
    alg_key << fZExpPath << "QEL-Gen0";
    r.Set(alg_key.str(), fZExpPara.fGen0);
    alg_key.str(""); // algorithm key for each coefficient
    alg_key << fZExpPath << "QEL-Gmn0";
    r.Set(alg_key.str(), fZExpPara.fGmn0);
    for (int i=0;i<fZExpParaDef.fKmax;i++)
    {
      alg_key.str(""); // algorithm key for each coefficient
      alg_key << fZExpPath << "QEL-Z_AN" << i+1;
      r.Set(alg_key.str(), fZExpPara.fZ_ANn[i]);
      alg_key.str(""); // algorithm key for each coefficient
      alg_key << fZExpPath << "QEL-Z_AP" << i+1;
      r.Set(alg_key.str(), fZExpPara.fZ_APn[i]);
      alg_key.str(""); // algorithm key for each coefficient
      alg_key << fZExpPath << "QEL-Z_BN" << i+1;
      r.Set(alg_key.str(), fZExpPara.fZ_BNn[i]);
      alg_key.str(""); // algorithm key for each coefficient
      alg_key << fZExpPath << "QEL-Z_BP" << i+1;
      r.Set(alg_key.str(), fZExpPara.fZ_BPn[i]);
    }
  }
  fXSecModel->Configure(r);
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEELFF::CalcWeight(const genie::EventRecord & event)
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

  double wght = 1.0;
      if ( fMode==kModeZExp && fModelIsZExp ) {
        wght *=  this->CalcWeightZExp( event );
        return wght;
      }

  return 1.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEELFF::Init(void)
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
  //  fFFModel = fXSecModelConfig->GetAlg("FormFactorsAlg/AxialFormFactorModel").name;
  fFFModel = fXSecModelConfig->GetAlg("FormFactorsAlg/ElasticFormFactorsModel").name;
  fXSecModelConfig->Print(std::cout);

  fModelIsZExp      = (strcmp(fFFModel.c_str(),kModelZExp  ) == 0);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  //this->SetMaPath("FormFactorsAlg/Ma");
  this->SetZExpPath("FormFactorsAlg/ElasticFormFactorsModel/");
  // this->SetZExpPath("FormFactorsAlg/AxialFormFactorModel/");

  if (fModelIsZExp)
  {
    this->SetMode(kModeZExp);
    fZExpParaDef.fQ4limit = fXSecModelConfig->GetBool(fZExpPath + "QEL-Q4limit");
    fZExpParaDef.fKmax    = fXSecModelConfig->GetInt(fZExpPath + "QEL-Kmax");
    fZExpParaDef.fT0      = fXSecModelConfig->GetDouble(fZExpPath + "QEL-T0");
    fZExpParaDef.fTcut    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Tcut");
    fZExpParaDef.fGep0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gep0");
    fZExpParaDef.fGmp0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gmp0");
    fZExpParaDef.fGen0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gen0");
    fZExpParaDef.fGmn0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gmn0");
    fZExpParaDef.fZ_ANn = new double [fZExpParaDef.fKmax];
    fZExpParaDef.fZ_APn = new double [fZExpParaDef.fKmax];
    fZExpParaDef.fZ_BNn = new double [fZExpParaDef.fKmax];
    fZExpParaDef.fZ_BPn = new double [fZExpParaDef.fKmax];
    fZExpPara.fZ_ANn = new double [fZExpParaDef.fKmax];
    fZExpPara.fZ_APn = new double [fZExpParaDef.fKmax];
    fZExpPara.fZ_BNn = new double [fZExpParaDef.fKmax];
    fZExpPara.fZ_BPn = new double [fZExpParaDef.fKmax];
    fZExpParaTwkDial.fZ_ANn = new double [fZExpParaDef.fKmax];
    fZExpParaTwkDial.fZ_APn = new double [fZExpParaDef.fKmax];
    fZExpParaTwkDial.fZ_BNn = new double [fZExpParaDef.fKmax];
    fZExpParaTwkDial.fZ_BPn = new double [fZExpParaDef.fKmax];
    fZExpParaTwkDial.fQ4limit = 0.;
    fZExpParaTwkDial.fKmax    = 0.;
    fZExpParaTwkDial.fT0      = 0.;
    fZExpParaTwkDial.fTcut    = 0.;
    fZExpParaTwkDial.fGep0    = 0.;
    fZExpParaTwkDial.fGmp0    = 0.;
    fZExpParaTwkDial.fGen0    = 0.;
    fZExpParaTwkDial.fGmn0    = 0.;



  ostringstream alg_key;
  for(int i = 0; i < fZExpParaDef.fKmax; i++){
    alg_key.str("");
    alg_key << fZExpPath << "QEL-Z_AN" << i+1;
    fZExpParaDef.fZ_ANn[i] = fXSecModelConfig->GetDouble(alg_key.str());
    alg_key.str("");
    alg_key << fZExpPath << "QEL-Z_AP" << i+1;
    fZExpParaDef.fZ_APn[i] = fXSecModelConfig->GetDouble(alg_key.str());
    alg_key.str("");
    alg_key << fZExpPath << "QEL-Z_BN" << i+1;
    fZExpParaDef.fZ_BNn[i] = fXSecModelConfig->GetDouble(alg_key.str());
    alg_key.str("");
    alg_key << fZExpPath << "QEL-Z_BP" << i+1;
    fZExpParaDef.fZ_BPn[i] = fXSecModelConfig->GetDouble(alg_key.str());
    fZExpParaTwkDial.fZ_ANn[i] = 0.;
    fZExpParaTwkDial.fZ_APn[i] = 0.;
    fZExpParaTwkDial.fZ_BNn[i] = 0.;
    fZExpParaTwkDial.fZ_BPn[i] = 0.;
    fZExpPara.fZ_ANn[i] = 0.;
    fZExpPara.fZ_APn[i] = 0.;
    fZExpPara.fZ_BNn[i] = 0.;
    fZExpPara.fZ_BPn[i] = 0.;
  }
  }

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  fTestFile = new TFile("./ccqe_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght");
#endif
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEELFF::CalcWeightZExp(const genie::EventRecord & event)
{
  // very similar to CalcWeightMa
  bool tweaked = false;
  tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fT0)   > controls::kASmallNum);
  tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fTcut) > controls::kASmallNum);
  tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fGep0) > controls::kASmallNum);
  tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fGmp0) > controls::kASmallNum);
  tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fGen0) > controls::kASmallNum);
  tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fGmn0) > controls::kASmallNum);
  for (int i=0;i<fZExpParaDef.fKmax;i++)
  {
    tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fZ_APn[i]) > controls::kASmallNum);
    tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fZ_ANn[i]) > controls::kASmallNum);
    tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fZ_BPn[i]) > controls::kASmallNum);
    tweaked = tweaked || (TMath::Abs(fZExpParaTwkDial.fZ_BNn[i]) > controls::kASmallNum);
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

