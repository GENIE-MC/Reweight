#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>

#include <TArrayF.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSpline.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

// GENIE/Generator includes
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

// GENIE/Reweight includes

#include "Rtypes.h"
#include "RwCalculators/GReWeightDeltaradAngle.h"
#include "RwCalculators/GReWeightINukeParams.h"
#include "RwCalculators/GReWeightNuXSecNC.h"
#include "RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "RwCalculators/GReWeightXSecMEC.h"

#include "RwCalculators/GReWeightProfessor.h"
#include "RwFramework/GReWeight.h"
#include "TArrayD.h"
#include "TAttLine.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TMatrixDfwd.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TVectorDfwd.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

#include <unordered_map>
#include <vector>

using namespace genie;
using namespace genie::rew;
using std::stringstream;

void PrintSyntax();
void GetEventRange(Long64_t nev_in_file, Long64_t &nfirst, Long64_t &nlast);
void GetCommandLineArgs(int argc, char **argv);
void GetCorrelationMatrix(string fname, TMatrixD *&cmat);
void AdoptWeightCalcs(vector<GSyst_t> lsyst, GReWeight &rw);
bool FindIncompatibleSystematics(vector<GSyst_t> lsyst);

// vector<GSyst_t> gOptVSyst; // not used as the syst are specifed in scanning
// file
static string gOptInpFilename;
static string gOptOutFilename;
static string gOptInpCovariance;
static Long64_t gOptNEvt1;
static Long64_t gOptNEvt2;
static int gOptRunKey = 0;
static int gOptNSyst = 0;
static int gOptNTwk = 0;
static long int gOptRanSeed; ///< random number seed

static std::string binningxml, ipol_path;

void GetCommandLineArgs(int argc, char **argv) {
  LOG("grwghtnp", pINFO) << "*** Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->ReadFromCommandLine(argc, argv);

  CmdLnArgParser parser(argc, argv);

  // get GENIE event sample
  if (parser.OptionExists('f')) {
    LOG("grwghtnp", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("grwghtnp", pFATAL) << "Unspecified input filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  // output weight file
  if (parser.OptionExists('o')) {
    LOG("grwghtnp", pINFO) << "Reading requested output filename";
    gOptOutFilename = parser.ArgAsString('o');
  } else {
    LOG("grwghtnp", pINFO) << "Setting default output filename";
    gOptOutFilename = "rwt_cov.root";
  }

  // get ROOT covariance matrix binary file
  if (parser.OptionExists('c')) {
    LOG("grwghtnp", pINFO) << "Reading covariance matrix filename";
    gOptInpCovariance = parser.ArgAsString('c');
  } else {
    LOG("grwghtnp", pFATAL) << "Unspecified covariance filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  if (parser.OptionExists('n')) {
    //
    LOG("grwghtnp", pINFO) << "Reading number of events to analyze";
    string nev = parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n', ",");
      if (vecn.size() != 2) {
        LOG("grwghtnp", pFATAL) << "Invalid syntax";
        gAbortingInErr = true;
        PrintSyntax();
        exit(1);
      }
      // User specified a comma-separated set of values n1,n2.
      // Use [n1,n2] as the event range to process.
      gOptNEvt1 = vecn[0];
      gOptNEvt2 = vecn[1];
    } else {
      // User specified a single number n.
      // Use [0,n] as the event range to process.
      gOptNEvt1 = -1;
      gOptNEvt2 = parser.ArgAsLong('n');
    }
  } else {
    LOG("grwghtnp", pINFO)
        << "Unspecified number of events to analyze - Use all";
    gOptNEvt1 = -1;
    gOptNEvt2 = -1;
  }
  LOG("grwghtnp", pDEBUG) << "Input event range: " << gOptNEvt1 << ", "
                          << gOptNEvt2;

  // number of tweaks:
  if (parser.OptionExists('t')) {
    LOG("grwghtnp", pINFO) << "Reading number of tweaks";
    gOptNTwk = parser.ArgAsInt('t');

    if (gOptNTwk < 1) {
      LOG("grwghtnp", pFATAL) << "Must have at least 1 tweak - Exiting";
      PrintSyntax();
      exit(1);
    }
    LOG("grwghtnp", pINFO) << "Number of tweaks : " << gOptNTwk;
  } else {
    LOG("grwghtnp", pFATAL) << "Unspecified tweaks for parameters - Exiting";
    PrintSyntax();
    exit(1);
  }

  // run key:
  if (parser.OptionExists('r')) {
    LOG("grwghtnp", pINFO) << "Reading run key";
    gOptRunKey = parser.ArgAsInt('r');

    LOG("grwghtnp", pINFO) << "Run key set to " << gOptRunKey;
  }
}

TMatrixD cov_to_corr(const TMatrixD &cov) {
  // Convert covariance matrix to correlation matrix
  TMatrixD corr(cov);
  for (int i = 0; i < cov.GetNrows(); ++i) {
    for (int j = 0; j < cov.GetNcols(); ++j) {
      if (i == j) {
        corr(i, j) = 1.0;
      } else {
        corr(i, j) = cov(i, j) / std::sqrt(cov(i, i) * cov(j, j));
      }
    }
  }
  return corr;
}

std::vector<double> cov_to_err(const TMatrixD &cov) {
  // Convert covariance matrix to a vector of uncertainties
  std::vector<double> err(cov.GetNrows());
  for (int i = 0; i < cov.GetNrows(); ++i) {
    err[i] = std::sqrt(cov(i, i));
  }
  return err;
}

int main(int argc, char **argv) {
  GetCommandLineArgs(argc, argv);
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  if (!RunOpt::Instance()->Tune()) {
    LOG("greweight", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // open the ROOT file and get the TTree & its header
  TTree *tree = 0;
  NtpMCTreeHeader *thdr = 0;
  TFile file(gOptInpFilename.c_str(), "READ");
  tree = dynamic_cast<TTree *>(file.Get("gtree"));
  thdr = dynamic_cast<NtpMCTreeHeader *>(file.Get("header"));
  if (!tree) {
    LOG("grwghtnp", pFATAL)
        << "Can't find a GHEP tree in input file: " << file.GetName();
    gAbortingInErr = true;
    PrintSyntax();
    exit(1);
  }
  LOG("grwghtnp", pNOTICE) << "Input tree header: " << *thdr;

  // LOG("grwghtnp", pNOTICE) << "Correlation matrix:";
  // cmat->Print();
  // LOG("grwghtnp", pNOTICE) << "Lower triangular matrix:";
  // lTri.Print();

  NtpMCEventRecord *mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  Long64_t nev_in_file = tree->GetEntries();
  Long64_t nfirst = 0;
  Long64_t nlast = 0;
  GetEventRange(nev_in_file, nfirst, nlast);
  int nev = int(nlast - nfirst + 1);

  auto observable_splines = std::make_unique<GReWeightProfessor>("");
  observable_splines->ReadComparisonXML(binningxml, ipol_path);

  std::vector<std::vector<double>> systematics_values{};

  // Gets Cor, which is needed in decompositions
  // Assumed errors from covariance are stored in one sigma errors for
  // parameters
  // GetCorrelationMatrix(gOptInpCovariance, cmat);
  TFile input_file = TFile(gOptInpCovariance.c_str(), "READ");
  auto cov_matrix = input_file.Get<TMatrixD>("param_covariance");
  if (!cov_matrix) {
    LOG("grwghtnp", pFATAL)
        << "Cannot find covariance matrix in file " << gOptInpCovariance;
    exit(1);
  }
  auto corr_matrix = cov_to_corr(*cov_matrix);
  auto std_err = cov_to_err(*cov_matrix);

  auto h_param_results =
      input_file.Get<TArrayD>("h_param_result"); // central values
  if (!h_param_results) {
    LOG("grwghtnp", pFATAL)
        << "Cannot find h_param_result in file " << gOptInpCovariance;
    exit(1);
  }

  std::vector<double> central_values{h_param_results->GetArray(),
                                     h_param_results->GetArray() +
                                         h_param_results->GetSize()};

  using namespace genie::utils::math;
  TMatrixD lTri = CholeskyDecomposition(corr_matrix);
  // todo: read the values from tunning file

  //

  for (int itk = 0; itk < gOptNTwk; itk++) {
    auto twkvals = CholeskyGenerateCorrelatedParamVariations(lTri);
    assert(twkvals.GetNrows() == central_values.size());
    auto vec =
        std::vector<double>(twkvals.GetMatrixArray(),
                            twkvals.GetMatrixArray() + twkvals.GetNrows());
    for (int i = 0; i < vec.size(); ++i) {
      // Apply the tweak to the central value
      vec[i] = central_values[i] + vec[i] * std_err[i];
    }
    systematics_values.push_back(vec);
  }

  auto *wght_file = new TFile(gOptOutFilename.c_str(), "RECREATE");
  auto *wght_tree = new TTree("covrwt", "GENIE covariant reweighting tree");
  wght_tree->SetDirectory(wght_file);
  TArrayD branch_weight(gOptNTwk);
  wght_tree->Branch("weights", &branch_weight);

  for (auto id_event{nfirst}; id_event <= nlast; ++id_event) {
    // Process each event with the corresponding systematic variations
    tree->GetEntry(id_event);
    EventRecord &event = *(mcrec->event);
    for (size_t idx = 0; idx < systematics_values.size(); ++idx) {
      const auto &sys_var = systematics_values[idx];
      // Apply each systematic variation to the event
      observable_splines->SetSystematic(sys_var, central_values);
      auto weight = observable_splines->CalcWeight(event);
      branch_weight[idx] = weight;
    }
    wght_tree->Fill();
    mcrec->Clear();
  }

  wght_file->Write();
  wght_file->Close();
  delete wght_file;

  return 0;
}