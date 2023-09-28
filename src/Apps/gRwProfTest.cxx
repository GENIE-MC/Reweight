#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <string>

#include <TArrayF.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
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
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

// GENIE/Reweight includes
#include "ROOT/RResultPtr.hxx"
#include "RwCalculators/GReWeightAGKY.h"
#include "RwCalculators/GReWeightDISNuclMod.h"
#include "RwCalculators/GReWeightFGM.h"
#include "RwCalculators/GReWeightFZone.h"
#include "RwCalculators/GReWeightINuke.h"
#include "RwCalculators/GReWeightNonResonanceBkg.h"
#include "RwCalculators/GReWeightNuXSecCCQE.h"
#include "RwCalculators/GReWeightNuXSecCCQEaxial.h"
#include "RwCalculators/GReWeightNuXSecCCQEvec.h"
#include "RwCalculators/GReWeightNuXSecCCRES.h"
#include "RwCalculators/GReWeightNuXSecCOH.h"
#include "RwCalculators/GReWeightNuXSecDIS.h"
#include "RwCalculators/GReWeightNuXSecNCEL.h"
#include "RwCalculators/GReWeightNuXSecNCRES.h"
#include "RwCalculators/GReWeightResonanceDecay.h"
#include "RwFramework/GReWeight.h"
#include "RwFramework/GReWeightI.h"
#include "RwFramework/GSyst.h"
#include "RwFramework/GSystSet.h"

#include "RwCalculators/GReWeightDeltaradAngle.h"
#include "RwCalculators/GReWeightINukeParams.h"
#include "RwCalculators/GReWeightNuXSecNC.h"
#include "RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "RwCalculators/GReWeightXSecMEC.h"

#include "ProfSpline/ObservableBins.h"
#include "RwCalculators/GReWeightProfessor.h"
#include "TAttLine.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSpline.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <unordered_map>
#include <vector>

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::rew;

std::vector<double> fromdat(std::string path) {
  std::vector<double> data{};
  std::ifstream infile(path);
  for (std::string line; std::getline(infile, line);) {
    std::istringstream iss(line);
    std::string _{};
    double a;
    iss >> _ >> a;
    data.push_back(a);
  }
  return data;
}

template <typename T>
std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
  TFile root_file{file_path.c_str(), "READ"};
  auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
  assert(objptr);
  return std::unique_ptr<T>{objptr};
}

std::pair<double, double> get_xsec(TH1 *h_rate, TGraph *spline) {
  double fluxint{};
  // spline->SaveAs("wrong.root");
  TSpline3 sp("sp", spline);
  TF1 func(
      "spline", [&](double *x, double *) { return sp.Eval(*x); }, 0,
      h_rate->GetXaxis()->GetXmax(), 0);
  for (int ii = 1; ii <= h_rate->GetNbinsX(); ii++) {
    double bin_c = h_rate->GetBinContent(ii);
    double bin_up = h_rate->GetXaxis()->GetBinUpEdge(ii);
    double bin_low = h_rate->GetXaxis()->GetBinLowEdge(ii);
    double bin_width = bin_up - bin_low;
    if (bin_c < 1 || func.Integral(bin_low, bin_up) == 0) {
      continue;
    }
    fluxint += bin_c / func.Integral(bin_low, bin_up) * bin_width;
  }
  double event_rate = h_rate->Integral();
  return {event_rate, event_rate / fluxint};
}

ROOT::RDF::RResultPtr<TH1D> normalize_filter(ROOT::RDF::RNode df) {
  return df
      .Filter(
          [](const NtpMCEventRecord &event) {
            return event.event->Summary()->ProcInfo().IsWeakCC() &&
                   event.event->Summary()->InitState().TgtPdg() == 1000060120;
          },
          {"gmcrec"})
      .Define("neutrinoE",
              [](const NtpMCEventRecord &event) {
                return event.event->Probe()->E();
              },
              {"gmcrec"})
      .Histo1D({"", "", 256, 0, 0}, "neutrinoE");
}

double get_genie_normalize(ROOT::RDF::RNode df, std::string filename, int Z) {
  auto h =
      df.Filter(
            [](const NtpMCEventRecord &event) {
              return event.event->Summary()->ProcInfo().IsWeakCC() &&
                     event.event->Summary()->InitState().TgtPdg() == 1000060120;
            },
            {"gmcrec"})
          .Define("neutrinoE",
                  [](const NtpMCEventRecord &event) {
                    return event.event->Probe()->E();
                  },
                  {"gmcrec"})
          .Histo1D({"", "", 256, 0, 0}, "neutrinoE");
  auto spline_obj = get_object<TGraph>(filename, "nu_mu_bar_H1/tot_cc");
  auto [tot, xsec] = get_xsec(h.GetPtr(), spline_obj.get());
  std::cerr << "tot: " << tot << " xsec: " << xsec << std::endl;
  xsec *= 1. / ((double)Z) * 1e-38;
  return xsec / tot;
}

template <typename T> auto common_def(T &&df) {
  return df
      .Filter(
          [](const NtpMCEventRecord &event) {
            // auto nucl = event.event->TargetNucleus();
            // if (!nucl)
            //   return false;
            return event.event->Summary()->ProcInfo().IsWeakCC() &&
                   event.event->Summary()->InitState().TgtPdg() == 1000060120;
          },
          {"gmcrec"})
      .Define("boost",
              [](const NtpMCEventRecord &event) {
                GHepParticle *target_nucleus_p = event.event->HitNucleon();

                auto target_nucleus = target_nucleus_p
                                          ? *target_nucleus_p
                                          : *(event.event->TargetNucleus());
                return -target_nucleus.P4()->BoostVector();
              },
              {"gmcrec"})
      .Define("muonp",
              [](const NtpMCEventRecord &event) {
                return event.event->FinalStatePrimaryLepton()->P4()->P();
              },
              {"gmcrec"})
      .Filter([](double muonp) { return muonp < 10.; }, {"muonp"})
      .Define("muon_theta",
              [](const NtpMCEventRecord &event) {
                return event.event->FinalStatePrimaryLepton()->P4()->Theta();
              },
              {"gmcrec"})
      .Define("q",
              [](const NtpMCEventRecord &event, const TVector3 &boost) {
                // return event.event->Summary()->KinePtr()->Q2();
                auto &p4_lep = *(event.event->FinalStatePrimaryLepton()->P4());
                auto &nu_lep = *(event.event->Probe()->P4());
                // GHepParticle *target_nucleus_p = event.event->HitNucleon();

                // auto target_nucleus = target_nucleus_p
                //                           ? *target_nucleus_p
                //                           : *(event.event->TargetNucleus());
                // auto boost_vec = -target_nucleus.P4()->BoostVector();
                auto vec = nu_lep - p4_lep;
                vec.Boost(boost);
                return vec;
              },
              {"gmcrec", "boost"})
      .Define("q0",
              [](const TLorentzVector &q) {
                // return event.event->Summary()->KinePtr()->Q2();
                return q.T();
              },
              {"q"})
      .Define("Q2",
              [](const TLorentzVector &q) {
                // return event.event->Summary()->KinePtr()->Q2();
                return -q.M2();
              },
              {"q"});
  // for now we filter out event outside of the binning range
}

// we don't have a proper configuration mechanism yet, so we need to
// hardcode sth.
int main(int argc, char **argv) {
  // std::vector<double> bin_edges_pmu{
  //     0,    0.5,  1,    1.5,  2,    2.5,  3,    3.5,  4,    4.5,  5,
  //     5.5,  6,    6.5,  7,    7.5,  8,    8.5,  9,    9.5,  10,   10.5,
  //     11,   11.5, 12,   12.5, 13,   13.5, 14,   14.5, 15,   15.5, 16,
  //     16.5, 17,   17.5, 18,   18.5, 19,   19.5, 20,   20.5, 100};
  // std::vector<double> bin_edges_enu{
  //     0,    0.5,  1,    1.5,  2,    2.5,  3,    3.5,  4,    4.5,  5,
  //     5.5,  6,    6.5,  7,    7.5,  8,    8.5,  9,    9.5,  10,   10.5,
  //     11,   11.5, 12,   12.5, 13,   13.5, 14,   14.5, 15,   15.5, 16,
  //     16.5, 17,   17.5, 18,   18.5, 19,   19.5, 20,   20.5, 100};
  // ObservableBins bins;
  // bins.InitializeBins({bin_edges_pmu, bin_edges_enu});
  // LOG("Test", pNOTICE) << bins.GetObservablesBinIDLinearized({2.1, 1.1});

  // ROOT::EnableImplicitMT();
  TH1::AddDirectory(false);
  std::string from_id = argc > 1 ? argv[1] : "0000";
  std::string to_id = argc > 2 ? argv[2] : "0002";
  std::string basedir = argc > 3 ? argv[3] : "/var/home/yan/neutrino/rewq2new/";
  if (basedir.back() != '/')
    basedir += "/";

  auto gevgen_filename = "minerva_antinumu_ccqe_NuMi.ghep.root";
  // auto gevgen_filename = "gntp.70000.ghep.root";

  std::string input_from_file = basedir + "scan/" + from_id +
                                "/master-professor_tune_01-minerva/" +
                                gevgen_filename;
  std::string input_to_file = basedir + "scan/" + to_id +
                              "/master-professor_tune_01-minerva/" +
                              gevgen_filename;
  auto spline_path_from = basedir + "scan/" + from_id +
                          "/master-professor_tune_01-xsec_vA/total_xsec.root";
  auto spline_path_to = basedir + "scan/" + to_id +
                        "/master-professor_tune_01-xsec_vA/total_xsec.root";
  auto ipol_path = basedir + "ipol_test.dat";

  ROOT::RDataFrame input_from_tree(
      "gtree",
      input_from_file); // read the tree from the file
  auto observable_splines = std::make_unique<GReWeightProfessor>("");
  observable_splines->Initialize("/var/home/yan/code/Comparisons/data/"
                                 "observables/pmu_enu_W/configure.xml");
  // std::string observable_name = "genie::rew::ObservablePMuEnu";
  // observable_splines->InitializeObservable(observable_name);
  // LOG("Test", pNOTICE) << "Initialize Professor Reweight success";

  observable_splines->ReadProf2Spline(ipol_path);
  LOG("Test", pNOTICE) << "Read Professor Reweight success";

  // // observable_splines->InitializeBins({bin_edges_pmu, bin_edges_enu});
  // LOG("Test", pNOTICE) << "Initialize Professor Reweight success";
  // now observable_splines should just work

  auto from = fromdat(basedir + "scan/" + from_id + "/params.dat");
  auto to = fromdat(basedir + "scan/" + to_id + "/params.dat");

  observable_splines->SetSystematic(to, from);
  LOG("Test", pNOTICE) << "Set Systematic success";

  LOG("Test", pNOTICE) << "Testing RDataFrame::Define to calculate weight\n";
  ROOT::RDataFrame tree_to{"gtree", input_to_file};

  auto ds_ref = common_def(tree_to);

  auto tree_rew = common_def(input_from_tree)
                      .Define("weight",
                              [&](const NtpMCEventRecord &event) {
                                auto weight = observable_splines->CalcWeight(
                                    *(event.event));
                                return weight < 0 ? 0 : weight;
                              },
                              {"gmcrec"});

  // auto hist_rew_alt = tree_rew.Histo1D(model, "muonp", "weight_alt");

  std::array<std::string, 4> varnames{"muonp", "muon_theta", "Q2", "q0"};

  std::unordered_map<std::string, ROOT::RDF::TH1DModel> modelmap{
      {"muonp", {"", ";p_{#mu};", 10, 0, 10.}},
      {"muon_theta", {"", ";#theta_{#mu};", 20, 0., 2.}},
      {"Q2", {"", ";Q^{2};", 30, 0., 6.}},
      {"q0", {"", ";q^{0};", 30, 0., 6.}}};

  // ROOT::RDF::TH1DModel model{"", "", 10, 0, 10.};

  std::unordered_map<
      std::string,
      std::tuple<ROOT::RDF::RResultPtr<TH1D>,
                 std::map<std::string, ROOT::RDF::RResultPtr<TH1D>>>>
      hist_rew_m, hist_reference_m, hist_unrew_m;

  std::unordered_map<std::string, std::function<bool(NtpMCEventRecord &)>>
      channelcuts{
          {"QEL",
           [](NtpMCEventRecord &event) {
             return event.event->Summary()->ProcInfo().IsQuasiElastic() &&
                    event.event->Summary()->ProcInfo().IsWeakCC();
           }},
          {"RES",
           [](NtpMCEventRecord &event) {
             return event.event->Summary()->ProcInfo().IsResonant() &&
                    event.event->Summary()->ProcInfo().IsWeakCC();
           }},
          {"DIS",
           [](NtpMCEventRecord &event) {
             return event.event->Summary()->ProcInfo().IsDeepInelastic() &&
                    event.event->Summary()->ProcInfo().IsWeakCC();
           }},
          {"2p2h",
           [](NtpMCEventRecord &event) {
             return event.event->Summary()->ProcInfo().IsMEC() &&
                    event.event->Summary()->ProcInfo().IsWeakCC();
           }},
          {"0pi",
           [](NtpMCEventRecord &event) {
             auto nentries = event.event->GetEntries();
             for (int i = 0; i < nentries; i++) {
               auto part = event.event->Particle(i);
               if (part->Status() == genie::kIStStableFinalState &&
                   (part->Pdg() == 211 || part->Pdg() == -211 ||
                    part->Pdg() == 111)) {
                 return false;
               }
             }
             return true;
           }},
          {"1pi",
           [](NtpMCEventRecord &event) {
             auto nentries = event.event->GetEntries();
             int count = 0;
             for (int i = 0; i < nentries; i++) {
               auto part = event.event->Particle(i);
               if (part->Status() == genie::kIStStableFinalState &&
                   (part->Pdg() == 211 || part->Pdg() == -211 ||
                    part->Pdg() == 111)) {
                 count++;
               }
             }
             return count == 1;
           }},
          {"Mpi", [](NtpMCEventRecord &event) {
             auto nentries = event.event->GetEntries();
             int count = 0;
             for (int i = 0; i < nentries; i++) {
               auto part = event.event->Particle(i);
               if (part->Status() == genie::kIStStableFinalState &&
                   (part->Pdg() == 211 || part->Pdg() == -211 ||
                    part->Pdg() == 111)) {
                 count++;
               }
             }
             return count > 1;
           }}};

  for (auto &&var : varnames) {
    auto &model = modelmap[var];
    std::get<0>(hist_rew_m[var]) =
        tree_rew.Histo1D<double, double>(model, var, "weight");
    std::get<0>(hist_reference_m[var]) = ds_ref.Histo1D<double>(model, var);
    std::get<0>(hist_unrew_m[var]) = tree_rew.Histo1D<double>(model, var);
    for (auto &&[channel, cut] : channelcuts) {
      std::get<1>(hist_rew_m[var])[channel] =
          tree_rew.Filter(cut, {"gmcrec"})
              .Histo1D<double, double>(model, var, "weight");
      std::get<1>(hist_reference_m[var])[channel] =
          ds_ref.Filter(cut, {"gmcrec"}).Histo1D<double>(model, var);
      std::get<1>(hist_unrew_m[var])[channel] =
          tree_rew.Filter(cut, {"gmcrec"}).Histo1D<double>(model, var);
    }
  }

  auto normalize_to = get_genie_normalize(tree_to, spline_path_to, 1);
  auto normalize_from =
      get_genie_normalize(input_from_tree, spline_path_from, 1);
  // auto normalize_to = std::bind(get_genie_normalize,
  // normalize_filter(tree_to),
  //                               spline_path_to, 1);
  // auto normalize_from =
  //     std::bind(get_genie_normalize, normalize_filter(input_from_tree),
  //               spline_path_from, 1);

  auto do_plot = [&](std::string varname, bool do_normalize,
                     std::string cut = "") {
    // auto hist_rew = hist_rew_m[varname].GetPtr();
    // auto hist_reference = hist_reference_m[varname].GetPtr();
    // auto hist_unrew = hist_unrew_m[varname].GetPtr();
    TH1D *hist_rew, *hist_reference, *hist_unrew;
    if (cut.empty()) {
      hist_rew = std::get<0>(hist_rew_m[varname]).GetPtr();
      hist_reference = std::get<0>(hist_reference_m[varname]).GetPtr();
      hist_unrew = std::get<0>(hist_unrew_m[varname]).GetPtr();
    } else {
      hist_rew = std::get<1>(hist_rew_m[varname])[cut].GetPtr();
      hist_reference = std::get<1>(hist_reference_m[varname])[cut].GetPtr();
      hist_unrew = std::get<1>(hist_unrew_m[varname])[cut].GetPtr();
    }
    std::string ytitle = "Count (Count * Weight)";
    if (do_normalize) {
      hist_reference->Scale(normalize_to, "width");
      hist_rew->Scale(normalize_from, "width");
      hist_unrew->Scale(normalize_from, "width");
      ytitle = "d#it{#sigma}/d";
      ytitle += hist_rew->GetXaxis()->GetTitle();
      ytitle += "(cm^{2}/MeV/nucleon)";
    }

    auto chi2 = hist_reference->Chi2Test(hist_rew, "CHI2/NDF");

    auto c1 = std::make_unique<TCanvas>();
    c1->SetLeftMargin(0.15);
    c1->Draw();
    // c1->Divide(1, 2);
    // auto pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    // auto pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    auto pad1 = std::make_unique<TPad>("pad1", "", 0, 0.3, 1, 1);
    auto pad2 = std::make_unique<TPad>("pad2", "", 0, 0, 1, 0.3);
    // c1->cd(1);
    pad1->cd();
    auto max = std::max({hist_reference->GetMaximum(), hist_rew->GetMaximum(),
                         hist_unrew->GetMaximum()});
    gStyle->SetOptStat(0);
    auto legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);

    hist_rew->SetLineColor(kBlue);
    hist_rew->SetMaximum(max * 1.1);
    hist_rew->GetYaxis()->SetTitle(ytitle.c_str());
    hist_rew->Draw("hist E");
    legend->AddEntry(hist_rew, "reweighted", "l");
    hist_unrew->SetLineColor(kGreen);
    hist_unrew->SetMaximum(max * 1.1);
    // hist_unrew->SetLineStyle(kDotted);
    hist_unrew->GetYaxis()->SetTitle(ytitle.c_str());
    hist_unrew->Draw("hist E same");
    legend->AddEntry(hist_unrew, "unreweighted", "l");
    hist_reference->SetMaximum(max * 1.1);
    hist_reference->SetLineColor(kRed);
    // hist_reference->SetLineStyle(kDotted);
    hist_reference->GetYaxis()->SetTitle(ytitle.c_str());
    hist_reference->Draw("hist E SAME");
    legend->AddEntry(hist_reference, "reference", "l");
    legend->SetHeader(("chi2/ndf: " + std::to_string(chi2)).c_str());
    legend->Draw("SAME");
    // c1->cd(2);
    pad2->cd();
    auto hist_diff_rew_ref =
        std::unique_ptr<TH1D>{static_cast<TH1D *>(hist_rew->Clone())};
    hist_diff_rew_ref->Add(hist_reference, -1);
    auto hist_diff_unrew_ref =
        std::unique_ptr<TH1D>{static_cast<TH1D *>(hist_unrew->Clone())};
    hist_diff_unrew_ref->Add(hist_reference, -1);
    // (hist_rew - hist_reference) / (hist_unrew - hist_reference);
    hist_diff_rew_ref->Divide(hist_diff_unrew_ref.get());
    // hist_diff_rew_ref->SetMaximum();
    hist_diff_rew_ref->GetYaxis()->SetTitle("rew - ref / unrew - ref");
    hist_diff_rew_ref->Draw("hist");
    c1->cd();
    pad1->Draw();
    pad2->Draw();
    c1->Draw();
    std::string filename = varname + "_" +
                           (do_normalize ? "normalized" : "unnormalized") +
                           cut + ".pdf";
    c1->SaveAs(filename.c_str());
  };
  tree_rew.Snapshot("out", "out.root", {"weight", "Q2"});

  // do_plot("muonp", true);
  // do_plot("muonp", false);
  // do_plot("muon_theta", true);
  // do_plot("muon_theta", false);
  // do_plot("Q2", true);
  // do_plot("Q2", false);
  for (auto &&var : varnames) {
    // do_plot(var, false);
    do_plot(var, true);
    do_plot(var, true, "QEL");
    do_plot(var, true, "RES");
    do_plot(var, true, "DIS");
    do_plot(var, true, "2p2h");
    do_plot(var, true, "0pi");
    do_plot(var, true, "1pi");
    do_plot(var, true, "Mpi");
  }

  return 0;
}