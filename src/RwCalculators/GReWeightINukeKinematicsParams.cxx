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

using namespace genie;
using namespace genie::rew;
using namespace genie::constants;

//___________________________________________________________________________
GReWeightINukeKinematicsParams::GReWeightINukeKinematicsParams():
  fNPwgt("pn"),
  fPPwgt("pp") {}

//___________________________________________________________________________
void GReWeightINukeKinematicsParams::SetTargetA(int A) {(void) A;} // nothing, for now

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
}

//___________________________________________________________________________
void GReWeightINukeKinematicsParams::Reset() {
  fNPwgt.Reset();
  fPPwgt.Reset();
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
  // std::cerr << "UNIVERSE: " << fUniverse << " VALID: " << ((int)!!fRewt) << " dw: " << ((!fRewt) ? 0. : fRewt->Evaluate(tlab, costhcm)) << std::endl;
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
      LOG("RWINukeKin", pERROR)
          << "Error reading INTRANUKE/hN data from: " << filename;
      return;
  }
  if(cols<2) {
    LOG("RWINukeKin", pERROR)
      << "Error reading INTRANUKE/hN data from: " << filename;
    LOG("RWINukeKin", pERROR)
      << "Too few columns: " << cols;
    return;
  }

  LOG("RWINukeKin", pINFO)
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

     LOG("RWINukeKin", pDEBUG)
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

