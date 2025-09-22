#include <cassert>
#include <cstdlib>

#include <TSystem.h>

// GENIE/Generator includes
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightINukeKinematicsParams.h"
#include "RwCalculators/GReWeightUtils.h"
#include "RwFramework/GSystUncertainty.h"

#include <gsl/gsl_integration.h>

using namespace genie;
using namespace genie::rew;
using namespace genie::constants;

//___________________________________________________________________________
GReWeightINukeKinematicsParams::GReWeightINukeKinematicsParams():
  fNPwgt("pn"),
  fPPwgt("pp"),
  fFixPiPro(false),
  fBiasPiPro(0.) {}

//___________________________________________________________________________
void GReWeightINukeKinematicsParams::SetTargetA(int A) {(void) A;} // nothing, for now

double ThreeBodyLorentzWeight(double E, double m1, double m2, double m3, double m23) {
  double p1cm = sqrt((E*E - (m1 + m23)*(m1 + m23))*(E*E - (m1 - m23)*(m1 - m23))) / (2*E);
  double p2cm2 = sqrt((m23*m23 - (m2+m3)*(m2+m3))*(m23*m23 - (m2-m3)*(m2-m3))) / (2*m23);
  return p1cm*p2cm2;
}

double ThreeBodyBias(double t1p, double bias) {
  return exp(bias*t1p);
}

double ThreeBodyLorentzWeightBias(double E, double m1, double m2, double m3, double m23, double t1p, double bias) {
  return ThreeBodyLorentzWeight(E, m1, m2, m3, m23)*ThreeBodyBias(t1p, bias);
}

// param should point to a list of 4 doubles: E, m1, m2, m3
double ParamsThreeBodyLorentzWeight(double m23, void *param) {
  double *M = (double *)param;
  return ThreeBodyLorentzWeight(M[0], M[1], M[2], M[3], m23);
}

double AvgThreeBodyLorentzWeight(double E, double m1, double m2, double m3) {
  double M[4] = {E, m1, m2, m3};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &ParamsThreeBodyLorentzWeight;
  F.params = (void*)M;
  double result, error;
  gsl_integration_qags(&F, m2+m3, E-m1, 0, 1.0e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / (E - m1 - m2 - m3);
}

// param should point to a list of 7 doubles: E, m1, m2, m3, m23, mp, bias
double ParamsThreeBodyLorentzWeightBiasInner(double costh1p, void *param) {
  double *M = (double *)param;
  double E = M[0];
  double m1 = M[1];
  double m2 = M[2];
  double m3 = M[3];
  double m23 = M[4];
  double mp = M[5];
  double bias = M[6];

  // calculate the energy transfer from the projectile to particle 1
  double p1 = sqrt((E*E - (m1 + m23)*(m1 + m23))*(E*E - (m1 - m23)*(m1 - m23))) / (2*E);
  double e1 = sqrt(p1*p1 + m1*m1);
  // assume equal masses in collision
  double ep = E / 2;
  double pp = sqrt(ep*ep - mp*mp);
  double t1p = m1*m1 + mp*mp - 2*(ep*e1 - pp*p1*costh1p); // == (pp - p1)**2

  return ThreeBodyLorentzWeightBias(E, m1, m2, m3, m23, t1p, bias);
}

// param should point to a list of 6 doubles: E, m1, m2, m3, mp, bias
double ParamsThreeBodyLorentzWeightBias(double m23, void *param) {
  double *M = (double *)param;
  double E = M[0];
  double m1 = M[1];
  double m2 = M[2];
  double m3 = M[3];
  double mp = M[4];
  double bias = M[5];

  double Minner[7] = {E, m1, m2, m3, m23, mp, bias};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &ParamsThreeBodyLorentzWeightBiasInner;
  F.params = (void*)Minner;
  double result, error;
  gsl_integration_qags(&F, -1, 1, 0, 1.0e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / 2; // divide out integration range
}

double AvgThreeBodyLorentzWeightBias(double E, double m1, double m2, double m3, double mp, double bias) {
  double M[6] = {E, m1, m2, m3, mp, bias};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &ParamsThreeBodyLorentzWeightBias;
  F.params = (void*)M;
  double result, error;
  gsl_integration_qags(&F, m2+m3, E-m1, 0, 1.0e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / (E - m1 - m2 - m3);
}

// param should point to a list of 5 doubles: E, m1, m23, mp, bias
double ParamsThreeBodyBiasInner(double costh1p, void *param) {
  double *M = (double *)param;
  double E = M[0];
  double m1 = M[1];
  double m23 = M[2];
  double mp = M[3];
  double bias = M[4];

  // calculate the energy transfer from the projectile to particle 1
  double p1 = sqrt((E*E - (m1 + m23)*(m1 + m23))*(E*E - (m1 - m23)*(m1 - m23))) / (2*E);
  double e1 = sqrt(p1*p1 + m1*m1);
  // assume equal masses in collision
  double ep = E / 2;
  double pp = sqrt(ep*ep - mp*mp);
  double t1p = m1*m1 + mp*mp - 2*(ep*e1 - pp*p1*costh1p); // == (pp - p1)**2

  return ThreeBodyBias(t1p, bias);
}

// param should point to a list of 4 doubles: E, m1, mp, bias
double ParamsThreeBodyBias(double m23, void *param) {
  double *M = (double *)param;
  double E = M[0];
  double m1 = M[1];
  double mp = M[2];
  double bias = M[3];

  double Minner[5] = {E, m1, m23, mp, bias};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &ParamsThreeBodyBiasInner;
  F.params = (void*)Minner;
  double result, error;
  gsl_integration_qags(&F, -1, 1, 0, 1.0e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / 2; // divide out integration range
}

double AvgThreeBodyBias(double E, double m1, double m2, double m3, double mp, double bias) {
  double M[4] = {E, m1, mp, bias};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &ParamsThreeBodyBias;
  F.params = (void*)M;
  double result, error;
  gsl_integration_qags(&F, m2+m3, E-m1, 0, 1.0e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / (E - m1 - m2 - m3);
}

double GReWeightINukeKinematicsParams::CalcPionWeight(TLorentzVector p, TLorentzVector f1, TLorentzVector f2, TLorentzVector f3) {
  double weight = 1.;

  double E = (f1 + f2 + f3).M();
  if (fFixPiPro && fBiasPiPro != 0.) {
    // Momentum transfer from parent to child nulceon
    double t1p = (p - f1).M2();

    // mass of 2-3 system
    double m23 = (f2 + f3).M();

    double bias_fix_weight = ThreeBodyLorentzWeightBias(E, f1.M(), f2.M(), f3.M(), m23, t1p, fBiasPiPro); 
    double ave_weight = AvgThreeBodyLorentzWeightBias(E, f1.M(), f2.M(), f3.M(), p.M(), fBiasPiPro);
    weight = bias_fix_weight / ave_weight;
  }
  else if (fBiasPiPro != 0.) {
    // Momentum transfer from parent to child nulceon
    double t1p = (p - f1).M2();
    double bias_weight = ThreeBodyBias(t1p, fBiasPiPro);
    double ave_weight = AvgThreeBodyBias(E, f1.M(), f2.M(), f3.M(), p.M(), fBiasPiPro);
    weight = bias_weight / ave_weight;
  }
  else if (fFixPiPro) {
    // mass of 2-3 system
    double m23 = (f2 + f3).M();

    // fix angles
    TVector3 boost_to_COM = -(f1+f2+f3).BoostVector();
    TLorentzVector f1_COM(f1);
    f1_COM.Boost(boost_to_COM);
    TVector3 z_dir = (f1+f2+f3).Vect().Unit();
    double costh1 = f1_COM.Vect().Unit().Dot(z_dir);
    double sinth1 = sqrt(1 - costh1*costh1);

    TVector3 boost_to_COM2 = -(f2+f3).BoostVector();
    TLorentzVector f2_COM2(f2);
    f2_COM2.Boost(boost_to_COM2);
    TVector3 z_dir2 = (f2+f3).Vect().Unit();
    double costh2 = f2_COM2.Vect().Unit().Dot(z_dir2);
    double sinth2 = sqrt(1 - costh2*costh2);

    double fix_weight = ThreeBodyLorentzWeight(E, f1.M(), f2.M(), f3.M(), m23)*sinth1*sinth2;
    double ave_weight = AvgThreeBodyLorentzWeight(E, f1.M(), f2.M(), f3.M()) * (2./M_PI) * (2./M_PI);
    weight = fix_weight / ave_weight;
  }

  return weight;
}

double GReWeightINukeKinematicsParams::CalcWeight(int np_pp, float tlab, float costhcm) {
  if (np_pp == 0) 
    return fNPwgt.CalcWeight(tlab, costhcm);

  return fPPwgt.CalcWeight(tlab, costhcm);
}

//___________________________________________________________________________
void GReWeightINukeKinematicsParams::SetTwkDial(GSyst_t syst, double val)
{
  if (syst == kINukeKinematicsTwkDial_NP_N) {
    fNPwgt.SetUniverse((int)val);
  }
  else if (syst == kINukeKinematicsTwkDial_PP_N) {
    fPPwgt.SetUniverse((int)val);
  }
  else if (syst == kINukeKinematicsFixPiPro) {
    if (val != 0.) fFixPiPro = true;
  }
  else if (syst == kINukeKinematicsPiProBias) {
    if (val < 0) LOG("ReW", pERROR) << "It is unphysical to bias pion production with a negative 'B' value. Therefore, only positive sigma values are physical.";

    GSystUncertainty * fracerr = GSystUncertainty::Instance();
    fBiasPiPro = val*fracerr->OneSigmaErr(kINukeKinematicsPiProBias);
  }
  else if (syst == kINukeKinematicsPiProBiaswFix) {
    GSystUncertainty * fracerr = GSystUncertainty::Instance();
    fBiasPiPro = val*fracerr->OneSigmaErr(kINukeKinematicsPiProBiaswFix);
    fFixPiPro = true;
  }
}

//___________________________________________________________________________
void GReWeightINukeKinematicsParams::Reset() {
  fNPwgt.Reset();
  fPPwgt.Reset();
  fFixPiPro = false;
  fBiasPiPro = 0.;
}

void GReWeightINukeKinematicsParams::ReBounce::Reset() {
  fRewt.reset();
  fCache.clear();
}

//___________________________________________________________________________
void GReWeightINukeKinematicsParams::ReBounce::SetUniverse(int u) {
  fUniverse = u;

  if (u == 0) { // set to CV
    fRewt.reset();
    return;
  }

  LoadUniverse(u);
  fRewt = fCache[u];
}

double GReWeightINukeKinematicsParams::ReBounce::CalcWeight(float tlab, float costhcm) {
  // std::cout << "UNIVERSE: " << fUniverse << " VALID: " << ((int)!!fRewt) << " dw: " << ((!fRewt) ? 0. : fRewt->Evaluate(tlab, costhcm)) << std::endl;
  if (!fRewt) return 1;

  return 1 + fRewt->Evaluate(tlab, costhcm);
}

// Helper function for reading hN cross section file
//____________________________________________________________________________
void ReadhNFile(
  string filename, double ke, int npoints, int & curr_point,
  double * costh_array, double * wgt_array, int cols)
{
  // open
  std::ifstream hN_stream(filename.c_str(), std::ios::in);
  if(!hN_stream.good()) {
      LOG("ReW", pERROR)
          << "Error reading INTRANUKE/hN data from: " << filename;
      return;
  }
  if(cols<2) {
    LOG("ReW", pERROR)
      << "Error reading INTRANUKE/hN data from: " << filename;
    LOG("ReW", pERROR)
      << "Too few columns: " << cols;
    return;
  }

  LOG("ReW", pINFO)
     << "Reading INTRANUKE/hN data from: " << filename;

  // skip initial comments
  char cbuf[501];
  hN_stream.getline(cbuf,400);
  hN_stream.getline(cbuf,400);
  hN_stream.getline(cbuf,400);

  // read
  double angle = 0;
  double wgt  = 0;
  double trash = 0;

  for(int ip = 0; ip < npoints; ip++) {
     hN_stream >> angle >> wgt;

     for(int ic = 0; ic < (cols-2); ic++) {
       hN_stream >> trash;
     }

     LOG("ReW", pDEBUG)
       << "Adding data point: (KE = " << ke << " MeV, angle = "
       << angle << ", sigma = " << wgt << ")";
     costh_array[ip] = TMath::Cos(angle*kPi/180.);
     wgt_array [curr_point] = wgt;
     curr_point++;
  }
}

std::shared_ptr<BLI2DNonUnifGrid> LoadWgts(std::string fname) {
  string data_dir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
             string(gSystem->Getenv("GINUKEHADRONDATA")) :
             string(gSystem->Getenv("GENIE")) + string("/data/syst/diff_ang_variations/");

  // Hard-coded parameters for hadron data, taken from
  // Physics/HadronTransport/INukeHadroData2018.cxx
  const int hN_wgt_nfiles = 20;
  const int hN_wgt_points_per_file = 21;
  const int hN_wgt_npoints = hN_wgt_points_per_file * hN_wgt_nfiles;
  
  double hN_wgt_energies[hN_wgt_nfiles] = {
    50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
    550, 600, 650, 700, 750, 800, 850, 900, 950, 1000
  };

  double hN_wgt_costh [hN_wgt_points_per_file] = {};
  double hN_wgt_wgt  [hN_wgt_npoints] = {};
  
  int ipoint=0;
  
  for(int ifile = 0; ifile < hN_wgt_nfiles; ifile++) {
    // build filename
    std::ostringstream hN_datafile;
    double ke = hN_wgt_energies[ifile];
    hN_datafile << data_dir << "/" << fname << ke << ".txt";
    // read data
    ReadhNFile(
      hN_datafile.str(), ke, hN_wgt_points_per_file,
      ipoint, hN_wgt_costh, hN_wgt_wgt,2);
  }//loop over files

  return std::shared_ptr<BLI2DNonUnifGrid>(new BLI2DNonUnifGrid(hN_wgt_nfiles,hN_wgt_points_per_file,
    hN_wgt_energies,hN_wgt_costh,hN_wgt_wgt));

}

GReWeightINukeKinematicsParams::ReBounce::~ReBounce() {}

GReWeightINukeKinematicsParams::ReBounce::ReBounce(std::string n):
  fUniverse(0), fName(n) {}

void GReWeightINukeKinematicsParams::ReBounce::LoadUniverse(int u) {
  // Do nothing if we already have the wgt
  if (fCache.count(u) > 0) return; 

  std::string file_name = "diff_ang_univ" + std::to_string(u) + "/" + fName + "/" + fName;
  fCache[u] = LoadWgts(file_name);
}

