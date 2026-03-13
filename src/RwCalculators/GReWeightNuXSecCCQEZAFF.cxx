//____________________________________________________________________________
/*
   Copyright (c) 2003-2025, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org

   Liang Liu <liangliu \at fnal.gov>
   Fermi National Accelerator Laboratory

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
#include "GReWeightNuXSecCCQEZAFF.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"
#include "TRandom3.h"
#include "TVectorD.h"
#include "TDecompChol.h"

using namespace genie;
using namespace genie::rew;
using std::ostringstream;

static const char* kModelZExp          = "genie::ZExpAxialFormFactorModel";

const int GReWeightNuXSecCCQEZAFF::kModeZExp;

//_______________________________________________________________________________________
GReWeightNuXSecCCQEZAFF::GReWeightNuXSecCCQEZAFF() :
  GReWeightModel("CCQE"),
  fManualModelName(),
  fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEZAFF::GReWeightNuXSecCCQEZAFF(std::string model, std::string type) :
  GReWeightModel("CCQE"),
  fManualModelName(model),
  fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEZAFF::~GReWeightNuXSecCCQEZAFF()
{
  if ( fXSecModelConfig ) delete fXSecModelConfig;
  if ( fXSecModel ) delete fXSecModel;
  if ( fXSecModelDef ) delete fXSecModelDef;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEZAFF::IsHandled(GSyst_t syst) const
{
  // read form factor model and compare to mode
  bool handle;

  switch(syst) {
    case ( kXSecTwkDial_ZExpZAFF ) :
      if(fMode==kModeZExp && fModelIsZExp){
        handle = true;
      }else {
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
bool GReWeightNuXSecCCQEZAFF::AppliesTo(const EventRecord &event) const
{
  auto type = event.Summary()->ProcInfo().ScatteringTypeId();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if (type==kScQuasiElastic && is_cc) {
    return true;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEZAFF::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst))
  {
    LOG("ReW",pWARN) << "Systematic " << GSyst::AsString(syst) << " is not handled for algorithm "
      << fFFModel << " and mode " << fMode;
    return;
  }
  switch(syst) {
    case (kXSecTwkDial_ZExpZAFF):
      fZExpTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEZAFF::Reset(void)
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
    fZExpPara.fZ_An[i] = fZExpParaDef.fZ_An[i];
    fZExpParaTwkDial.fZ_An[i] = 0.;
  }
  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEZAFF::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();
  if(fMode==kModeZExp && fModelIsZExp) {
    int     sign_twk = 0;
    double  fracerr_zexp = 0.;
    sign_twk = utils::rew::Sign(fZExpTwkDial);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpZAFF, sign_twk);
    fZExp_Scale = fZExpTwkDial * fracerr_zexp;
  }
  else {
    return;
  }
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEZAFF::CalcWeight(const genie::EventRecord & event)
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
void GReWeightNuXSecCCQEZAFF::Init(void)
{

  // Get the model and parameters of axial form factor from current tune
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  // get axial form factor tune from current CCQE model
  RgAlg cc_qel_id = gpl->GetAlg( "XSecModel@genie::EventGenerator/QEL-CC" );
  RgAlg ff_id = conf_pool->FindRegistry(cc_qel_id)->GetAlg("FormFactorsAlg");
  RgAlg aff_model_id = conf_pool->FindRegistry(ff_id)->GetAlg("AxialFormFactorModel");
  // get the covariance matrix
  Registry * zexp_axial_model = conf_pool->FindRegistry(aff_model_id);
  int n_row = zexp_axial_model->GetInt(Algorithm::BuildParamMatRowSizeKey("ZExpZAFF@CovarianceMatrix"));
  int n_col = zexp_axial_model->GetInt(Algorithm::BuildParamMatColSizeKey("ZExpZAFF@CovarianceMatrix"));

  if(n_row != n_col){
    LOG( "GReWeightNuXSecCCQEZAFF", pFATAL ) << "Non-square covariance matrix"
      << "encountered in GReWeightNuXSecCCQEZAFF::Init()";
    std::exit(1);
  }
  error_mat.ResizeTo(n_row, n_row);
  errors.resize(n_row);
  A_f.resize(n_row);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_row; j++){
      error_mat[i][j] = zexp_axial_model->GetDouble(Algorithm::BuildParamMatKey("ZExpZAFF@CovarianceMatrix", i, j));
    }
  }
  LOG( "GReWeightNuXSecCCQEZAFF", pINFO ) << "CovarianceMatrix of z expansion:";
  error_mat.Print();

  AlgId id(cc_qel_id);

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


  // Check what kind of form factors we're using in the tweaked cross section
  // model
  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
  fFFModel = fXSecModelConfig->GetAlg("FormFactorsAlg/AxialFormFactorModel").name;
  fXSecModelConfig->Print(std::cout);

  fModelIsZExp      = (strcmp(fFFModel.c_str(), kModelZExp  ) == 0);


  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  this->SetZExpPath("FormFactorsAlg/AxialFormFactorModel/");

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

    fZExpParaDef.fZ_An.resize(fZExpParaDef.fKmax);
    fZExpPara.fZ_An.resize(fZExpParaDef.fKmax);
    fZExpParaTwkDial.fZ_An.resize(fZExpParaDef.fKmax);

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
      alg_key << fZExpPath << "QEL-Z_A-" << i;
      fZExpParaDef.fZ_An[i] = fXSecModelConfig->GetDouble(alg_key.str());
      fZExpParaTwkDial.fZ_An[i] = 0.;
      fZExpPara.fZ_An[i] = 0.;
    }
    fZExpTwkDial = 0.;
  }
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEZAFF::CalcWeightZExp(const genie::EventRecord & event)
{
    bool tweaked = false;
    tweaked = tweaked || (TMath::Abs(fZExpTwkDial) > controls::kASmallNum);
    if(!tweaked) { return 1.0; }
    double oneSigma = GetOneSigma(event);
    double old_xsec = event.DiffXSec();
    double new_weight = (fZExp_Scale * oneSigma + old_xsec) / old_xsec;
    return new_weight;
}

//_______________________________________________________________________________________
// Calculate the partial derivative of the cross section
// A_f is the partial derivative
// A_f = \patial XSec() / \patial p_i = (XSec(p_i + \delta * error_i ) - XSec(p_i - \delta * error_i )) / ( 2.0 * \delta * error_i )
//

void GReWeightNuXSecCCQEZAFF::XSecPartialDerivative(const EventRecord & event){
  // Get the uncertainties from the error matrix
  // ap1, ap2, ap3, ap4,
  // bp1, bp2, bp3, bp4,
  // an1, an2, an3, an4,
  // bn1, bn2, bn3, bn4

  for(int i = 0; i < fZExpParaDef.fKmax; i++){
    errors[i] = TMath::Sqrt(error_mat[i][i]);
    A_f[i] = 0.0;
  }

  double delta = 0.1;

  for(int index = 0; index < fZExpParaDef.fKmax; index++){
    double xsec_tmp_0 = 0.0;
    double xsec_tmp_1 = 0.0;

    for(int sign = -1; sign < 2; sign++){
      if(sign == 0) continue;

      for(int ipara = 0; ipara < fZExpParaDef.fKmax; ipara++){
        fZExpPara.fZ_An[ipara] = fZExpParaDef.fZ_An[ipara];
      }

      fZExpPara.fZ_An[index] = fZExpParaDef.fZ_An[index] + errors[index] * delta * sign; 

      Registry r("GReWeightNuXSecCCQEZAFF",false);
      //~ Registry r(fXSecModel->GetConfig());
      if (fMode==kModeZExp)
      {
        ostringstream alg_key;
        for (int i=0;i<fZExpParaDef.fKmax;i++)
        {
          alg_key.str(""); // algorithm key for each coefficient
          alg_key << fZExpPath << "QEL-Z_A-" << i;
          r.Set(alg_key.str(), fZExpPara.fZ_An[i]);
          LOG("ReW", pINFO) << alg_key.str() << "  " << fZExpParaDef.fZ_An[i] - fZExpPara.fZ_An[i] << "  " << fZExpPara.fZ_An[i];
        }
      }
      fXSecModel->Configure(r);
      if(sign == -1){
        xsec_tmp_0 = UpdateXSec(event) ;
      }
      else if(sign == 1){
        xsec_tmp_1 = UpdateXSec(event) ;
      }
    }

    double delta_xsec = xsec_tmp_1 - xsec_tmp_0;
    A_f[index] = delta_xsec/(errors[index] * delta * 2.0);
    //   LOG("ReW", pNOTICE) << A_f[index] ;
  }
}

//_______________________________________________________________________________________
//  Uncertainty propagation
//
//  \sigma_{XSec}^2 = A_f[i] *A_f[j] *M_ij
double GReWeightNuXSecCCQEZAFF::GetOneSigma(const EventRecord & event){
  XSecPartialDerivative(event);
  double OneSigma2 = 0;
  for(int i = 0; i < fZExpParaDef.fKmax; i++){
    for(int j = 0; j < fZExpParaDef.fKmax; j++){
      OneSigma2+= A_f[i]*A_f[j]*error_mat[i][j];
    }
  }
  return TMath::Sqrt(OneSigma2);
}

//_______________________________________________________________________________________


double GReWeightNuXSecCCQEZAFF::UpdateXSec(const EventRecord & event){

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
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);

  interaction->KinePtr()->ClearRunningValues();

  if (phase_space == kPSQ2fE) {
    interaction->ResetBit(kIAssumeFreeNucleon);
  }

  return new_xsec;
}


