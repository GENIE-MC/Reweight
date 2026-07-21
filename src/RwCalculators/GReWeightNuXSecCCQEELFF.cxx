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
#include "GReWeightNuXSecCCQEELFF.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"
#include "TRandom3.h"
#include "TVectorD.h"
#include "Framework/Numerical/MathUtils.h"
#include "TDecompChol.h"

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
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEELFF::IsHandled(GSyst_t syst) const
{
  // read form factor model and compare to mode
  bool handle;

  switch(syst) {
    // add ZExp vector parameters
    case ( kXSecTwkDial_ZExpELFF ) :
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
    case (kXSecTwkDial_ZExpELFF):
      fZExpTwkDial = twk_dial;
      break;
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
    fZExpParaTwkDial.fZ_APn[i] = 0.;
    fZExpParaTwkDial.fZ_BPn[i] = 0.;
    fZExpParaTwkDial.fZ_ANn[i] = 0.;
    fZExpParaTwkDial.fZ_BNn[i] = 0.;
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
    sign_twk = utils::rew::Sign(fZExpTwkDial);
    fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpELFF, sign_twk);
    fZExp_Scale = fZExpTwkDial * fracerr_zexp;
  }
  else {
    return;
  }
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

  // Finite-difference step for XSecPartialDerivative: reweight-tool numerical
  // setting (not model physics) — built-in default, per-job override via
  // SetFiniteDiffDelta() (exposed as `grwght1p --fd-delta`).
  fFiniteDiffDelta = 0.1;
  LOG( "GReWeightNuXSecCCQEELFF", pINFO )
    << "Finite-difference delta default: " << fFiniteDiffDelta
    << " (override with SetFiniteDiffDelta / grwght1p --fd-delta)";

  // sigma-estimator defaults (driver overrides via SetSigmaEstimator /
  // SetNUniverses, exposed as grwght1p --sigma-method / --n-universes)
  fSigmaEstimator = kSigmaPropagation;
  fNUniverses     = 1000;
  fLchComputed    = false;

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

  this->SetZExpPath("FormFactorsAlg/ElasticFormFactorsModel/");

  if (fModelIsZExp)
  {
    this->SetMode(kModeZExp);

    // Covariance matrix of the zexp elastic FF: only read it when the tune
    // actually uses the zexp EL model — reading unconditionally makes Init()
    // die (missing Registry key) on every tune without ZExpELFF@CovarianceMatrix,
    // which breaks ALL grwght1p dials, not just the ELFF one.
    RgAlg elff_id = conf_pool->FindRegistry("CommonParamList/ElasticFF")->GetAlg("ElasticFormFactorsModel");
    Registry * elff_model = conf_pool->FindRegistry(elff_id);
    int n_row = elff_model->GetInt(Algorithm::BuildParamMatRowSizeKey("ZExpELFF@CovarianceMatrix"));
    int n_col = elff_model->GetInt(Algorithm::BuildParamMatColSizeKey("ZExpELFF@CovarianceMatrix"));
    if(n_row != n_col){
      LOG( "GReWeightNuXSecCCQEELFF", pFATAL ) << "Non-square covariance matrix"
        << "encountered in GReWeightNuXSecCCQEELFF::Init()";
      std::exit(1);
    }
    error_mat.ResizeTo(n_row, n_row);
    errors.resize(n_row);
    A_f.resize(n_row);
    for(int i = 0; i < n_row; i++){
      for(int j = 0; j < n_row; j++){
        error_mat[i][j] = elff_model->GetDouble(Algorithm::BuildParamMatKey("ZExpELFF@CovarianceMatrix", i, j));
      }
    }
    error_mat.Print();

    fZExpParaDef.fQ4limit = fXSecModelConfig->GetBool(fZExpPath + "QEL-Q4limit");
    fZExpParaDef.fKmax    = fXSecModelConfig->GetInt(fZExpPath + "QEL-Kmax");
    fZExpParaDef.fT0      = fXSecModelConfig->GetDouble(fZExpPath + "QEL-T0");
    fZExpParaDef.fTcut    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Tcut");
    fZExpParaDef.fGep0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gep0");
    fZExpParaDef.fGmp0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gmp0");
    fZExpParaDef.fGen0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gen0");
    fZExpParaDef.fGmn0    = fXSecModelConfig->GetDouble(fZExpPath + "QEL-Gmn0");

    // Cross-checks (mirror of GReWeightNuXSecCCQEZAFF): the covariance is over
    // the 4 stacked coefficient blocks (AP,BP,AN,BN) x Kmax, and each QEL-Z_*
    // vector must have exactly Kmax entries — otherwise the derivative loops
    // in XSecPartialDerivative()/GetOneSigma() index out of bounds.
    if ( 4 * fZExpParaDef.fKmax != error_mat.GetNrows() ) {
      LOG( "GReWeightNuXSecCCQEELFF", pFATAL )
        << "4 x QEL-Kmax (" << 4 * fZExpParaDef.fKmax
        << ") does not match the ZExpELFF@CovarianceMatrix dimension ("
        << error_mat.GetNrows() << "x" << error_mat.GetNcols()
        << ") for elastic form factor model " << fFFModel
        << " - fix the model configuration (config set / GXMLPATH override)";
      std::exit(1);
    }
    const char * elff_vec_names[4] = { "QEL-Z_AP", "QEL-Z_BP", "QEL-Z_AN", "QEL-Z_BN" };
    for (int iv = 0; iv < 4; ++iv) {
      string vec_size_key = fZExpPath + Algorithm::BuildParamVectSizeKey(elff_vec_names[iv]);
      if ( ! fXSecModelConfig->Exists(RgKey(vec_size_key)) ) {
        LOG( "GReWeightNuXSecCCQEELFF", pFATAL )
          << "Missing '" << vec_size_key << "' - " << elff_vec_names[iv]
          << " must be configured as a vec-double parameter";
        std::exit(1);
      }
      int n_entries = fXSecModelConfig->GetInt(vec_size_key);
      if ( n_entries != fZExpParaDef.fKmax ) {
        LOG( "GReWeightNuXSecCCQEELFF", pFATAL )
          << elff_vec_names[iv] << " vector length (" << n_entries
          << ") does not match QEL-Kmax (" << fZExpParaDef.fKmax
          << ") for elastic form factor model " << fFFModel;
        std::exit(1);
      }
    }

    fZExpParaDef.fZ_ANn.resize(fZExpParaDef.fKmax);
    fZExpParaDef.fZ_APn.resize(fZExpParaDef.fKmax);
    fZExpParaDef.fZ_BNn.resize(fZExpParaDef.fKmax);
    fZExpParaDef.fZ_BPn.resize(fZExpParaDef.fKmax);
    fZExpPara.fZ_ANn.resize(fZExpParaDef.fKmax);
    fZExpPara.fZ_APn.resize(fZExpParaDef.fKmax);
    fZExpPara.fZ_BNn.resize(fZExpParaDef.fKmax);
    fZExpPara.fZ_BPn.resize(fZExpParaDef.fKmax);
    fZExpParaTwkDial.fZ_ANn.resize(fZExpParaDef.fKmax);
    fZExpParaTwkDial.fZ_APn.resize(fZExpParaDef.fKmax);
    fZExpParaTwkDial.fZ_BNn.resize(fZExpParaDef.fKmax);
    fZExpParaTwkDial.fZ_BPn.resize(fZExpParaDef.fKmax);

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
      alg_key << fZExpPath << "QEL-Z_AN-" << i;
      fZExpParaDef.fZ_ANn[i] = fXSecModelConfig->GetDouble(alg_key.str());
      alg_key.str("");
      alg_key << fZExpPath << "QEL-Z_AP-" << i;
      fZExpParaDef.fZ_APn[i] = fXSecModelConfig->GetDouble(alg_key.str());
      alg_key.str("");
      alg_key << fZExpPath << "QEL-Z_BN-" << i;
      fZExpParaDef.fZ_BNn[i] = fXSecModelConfig->GetDouble(alg_key.str());
      alg_key.str("");
      alg_key << fZExpPath << "QEL-Z_BP-" << i;
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
    fZExpTwkDial = 0.;
  }
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEELFF::CalcWeightZExp(const genie::EventRecord & event)
{
    bool tweaked = false;
    tweaked = tweaked || (TMath::Abs(fZExpTwkDial) > controls::kASmallNum);
    if(!tweaked) { return 1.0; }
    double oneSigma = GetOneSigma(event);
    double old_xsec = event.DiffXSec();
    double new_weight = (fZExp_Scale * oneSigma + old_xsec) / old_xsec;
 //   LOG("ReW", pNOTICE) << fZExp_Scale << "  " << oneSigma << "  " << old_xsec;
    return new_weight;
}

//_______________________________________________________________________________________
// Calculate the partial derivative of the cross section
// A_f is the partial derivative
// A_f = \patial XSec() / \patial p_i = (XSec(p_i + \delta * error_i ) - XSec(p_i - \delta * error_i )) / ( 2.0 * \delta * error_i )
//

void GReWeightNuXSecCCQEELFF::XSecPartialDerivative(const EventRecord & event){
  // Get the uncertainties from the error matrix
  // ap1, ap2, ap3, ap4,
  // bp1, bp2, bp3, bp4,
  // an1, an2, an3, an4,
  // bn1, bn2, bn3, bn4

  if ( fFiniteDiffDelta <= 0. ) {
    LOG( "GReWeightNuXSecCCQEELFF", pFATAL )
      << "Finite-difference delta must be > 0, got " << fFiniteDiffDelta;
    std::exit(1);
  }

  for(int i = 0; i < fZExpParaDef.fKmax * 4; i++){
    errors[i] = TMath::Sqrt(error_mat[i][i]);
    A_f[i] = 0.0;
  }

  double delta = fFiniteDiffDelta;

  for(int index = 0; index < fZExpParaDef.fKmax * 4; index++){
    double xsec_tmp_0 = 0.0;
    double xsec_tmp_1 = 0.0;

    for(int sign = -1; sign < 2; sign++){
      if(sign == 0) continue;

      for(int ipara = 0; ipara < fZExpParaDef.fKmax; ipara++){
        fZExpPara.fZ_APn[ipara] = fZExpParaDef.fZ_APn[ipara];
        fZExpPara.fZ_BPn[ipara] = fZExpParaDef.fZ_BPn[ipara];
        fZExpPara.fZ_ANn[ipara] = fZExpParaDef.fZ_ANn[ipara];
        fZExpPara.fZ_BNn[ipara] = fZExpParaDef.fZ_BNn[ipara];
      }

      int icoff = index / fZExpParaDef.fKmax;
      int jcoff = index % fZExpParaDef.fKmax;

      switch (icoff){
        case 0 : fZExpPara.fZ_APn[jcoff] = fZExpParaDef.fZ_APn[jcoff] + errors[index] * delta * sign; break;
        case 1 : fZExpPara.fZ_BPn[jcoff] = fZExpParaDef.fZ_BPn[jcoff] + errors[index] * delta * sign; break;
        case 2 : fZExpPara.fZ_ANn[jcoff] = fZExpParaDef.fZ_ANn[jcoff] + errors[index] * delta * sign; break;
        case 3 : fZExpPara.fZ_BNn[jcoff] = fZExpParaDef.fZ_BNn[jcoff] + errors[index] * delta * sign; break;
        default: break;
      }

      Registry r("GReWeightNuXSecCCQEELFF",false);
      //~ Registry r(fXSecModel->GetConfig());
      if (fMode==kModeZExp)
      {
        ostringstream alg_key;
        for (int i=0;i<fZExpParaDef.fKmax;i++)
        {
          alg_key.str(""); // algorithm key for each coefficient
          alg_key << fZExpPath << "QEL-Z_AN-" << i;
          r.Set(alg_key.str(), fZExpPara.fZ_ANn[i]);
          LOG("ReW", pINFO) << alg_key.str() << "  " << fZExpParaDef.fZ_ANn[i] - fZExpPara.fZ_ANn[i] << "  " << fZExpPara.fZ_ANn[i];
          alg_key.str(""); // algorithm key for each coefficient
          alg_key << fZExpPath << "QEL-Z_AP-" << i;
          r.Set(alg_key.str(), fZExpPara.fZ_APn[i]);
          LOG("ReW", pINFO) << alg_key.str() << "  " << fZExpParaDef.fZ_APn[i] - fZExpPara.fZ_APn[i] << "  " << fZExpPara.fZ_APn[i];
          alg_key.str(""); // algorithm key for each coefficient
          alg_key << fZExpPath << "QEL-Z_BN-" << i;
          r.Set(alg_key.str(), fZExpPara.fZ_BNn[i]);
          LOG("ReW", pINFO) << alg_key.str() << "  " << fZExpParaDef.fZ_BNn[i] - fZExpPara.fZ_BNn[i] << "  " << fZExpPara.fZ_BNn[i];
          alg_key.str(""); // algorithm key for each coefficient
          alg_key << fZExpPath << "QEL-Z_BP-" << i;
          r.Set(alg_key.str(), fZExpPara.fZ_BPn[i]);
          LOG("ReW", pINFO) << alg_key.str() << "  " << fZExpParaDef.fZ_BPn[i] - fZExpPara.fZ_BPn[i] <<  "  " <<  fZExpPara.fZ_BPn[i];
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
double GReWeightNuXSecCCQEELFF::GetOneSigma(const EventRecord & event){
  if ( fSigmaEstimator == kSigmaCholesky ) {
    return this->GetOneSigmaCholesky(event);
  }
  XSecPartialDerivative(event);
  double OneSigma2 = 0;
  for(int i = 0; i < fZExpParaDef.fKmax * 4; i++){
    for(int j = 0; j < fZExpParaDef.fKmax * 4; j++){
      OneSigma2+= A_f[i]*A_f[j]*error_mat[i][j];
    }
  }
  return TMath::Sqrt(OneSigma2);
}

//_______________________________________________________________________________________
//  Alternative sigma estimator: Cholesky-sampled universes (mirror of
//  GReWeightNuXSecCCQEZAFF::GetOneSigmaCholesky). The covariance — and hence
//  the variation vector — spans the 4*Kmax stacked coefficient blocks in the
//  order (AP, BP, AN, BN), matching XSecPartialDerivative's index mapping.
double GReWeightNuXSecCCQEELFF::GetOneSigmaCholesky(const EventRecord & event){
  if ( fNUniverses < 2 ) {
    LOG( "GReWeightNuXSecCCQEELFF", pFATAL )
      << "Number of universes must be >= 2, got " << fNUniverses;
    std::exit(1);
  }

  // Lazily compute (and validate) the Cholesky factor on first use, so the
  // propagation estimator never depends on the covariance being
  // positive-definite.
  if ( !fLchComputed ) {
    TDecompChol chol( error_mat );
    if ( !chol.Decompose() ) {
      LOG( "GReWeightNuXSecCCQEELFF", pFATAL )
        << "ZExpELFF@CovarianceMatrix is not positive-definite - cannot use "
        << "the Cholesky sigma estimator. Fix the covariance data or use "
        << "--sigma-method propagation.";
      std::exit(1);
    }
    fLch.ResizeTo(error_mat.GetNrows(), error_mat.GetNcols());
    fLch = TMatrixD(TMatrixD::kTransposed, chol.GetU());  // lower-triangular L
    fLchComputed = true;
  }

  const int kmax = fZExpParaDef.fKmax;
  double sum = 0., sum2 = 0.;
  ostringstream alg_key;
  for (int u = 0; u < fNUniverses; ++u) {
    TVectorD var = genie::utils::math::CholeskyGenerateCorrelatedParamVariations(fLch);
    // stacked ordering: index = icoff*kmax + jcoff, icoff 0=AP, 1=BP, 2=AN, 3=BN
    for (int j = 0; j < kmax; ++j) {
      fZExpPara.fZ_APn[j] = fZExpParaDef.fZ_APn[j] + var[0*kmax + j];
      fZExpPara.fZ_BPn[j] = fZExpParaDef.fZ_BPn[j] + var[1*kmax + j];
      fZExpPara.fZ_ANn[j] = fZExpParaDef.fZ_ANn[j] + var[2*kmax + j];
      fZExpPara.fZ_BNn[j] = fZExpParaDef.fZ_BNn[j] + var[3*kmax + j];
    }
    Registry r("GReWeightNuXSecCCQEELFF",false);
    for (int i = 0; i < kmax; ++i) {
      alg_key.str(""); alg_key << fZExpPath << "QEL-Z_AN-" << i; r.Set(alg_key.str(), fZExpPara.fZ_ANn[i]);
      alg_key.str(""); alg_key << fZExpPath << "QEL-Z_AP-" << i; r.Set(alg_key.str(), fZExpPara.fZ_APn[i]);
      alg_key.str(""); alg_key << fZExpPath << "QEL-Z_BN-" << i; r.Set(alg_key.str(), fZExpPara.fZ_BNn[i]);
      alg_key.str(""); alg_key << fZExpPath << "QEL-Z_BP-" << i; r.Set(alg_key.str(), fZExpPara.fZ_BPn[i]);
    }
    fXSecModel->Configure(r);
    double xs = this->UpdateXSec(event);
    sum  += xs;
    sum2 += xs*xs;
  }
  double mean = sum / fNUniverses;
  double var_x = (sum2 - fNUniverses*mean*mean) / (fNUniverses - 1.);
  if ( var_x < 0. ) var_x = 0.;  // numerical guard
  return TMath::Sqrt(var_x);
}

//_______________________________________________________________________________________


double GReWeightNuXSecCCQEELFF::UpdateXSec(const EventRecord & event){

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


