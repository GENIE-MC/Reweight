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
#include "Framework/Utils/XSecSplineList.h"

// GENIE/Reweight includes

#include "RwCalculators/GReWeightDeltaradAngle.h"
#include "RwCalculators/GReWeightINukeParams.h"
#include "RwCalculators/GReWeightNuXSecNC.h"
#include "RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "RwCalculators/GReWeightXSecMEC.h"

#include "RwCalculators/GReWeightProfessor.h"
#include "TAttLine.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSpline.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
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

template <typename T> auto common_def(T &&df) {
  return df
      .Filter(
          [](const NtpMCEventRecord &event) {
            return event.event->Summary()->ProcInfo().IsWeakCC() &&
                   (event.event->Summary()->ProcInfo().IsResonant() ||
                    event.event->Summary()->ProcInfo().IsDeepInelastic())

                ;
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
      .Define("muon_theta",
              [](const NtpMCEventRecord &event) {
                return event.event->FinalStatePrimaryLepton()->P4()->Theta();
              },
              {"gmcrec"})
      .Define("q",
              [](const NtpMCEventRecord &event, const TVector3 &boost) {
                // return event.event.event->Summary()->KinePtr()->Q2();
                auto &p4_lep = *(event.event->FinalStatePrimaryLepton()->P4());
                auto &nu_lep = *(event.event->Probe()->P4());
                // GHepParticle *target_nucleus_p =
                // event.event.event->HitNucleon();

                // auto target_nucleus = target_nucleus_p
                //                           ? *target_nucleus_p
                //                           :
                //                           *(event.event.event->TargetNucleus());
                // auto boost_vec = -target_nucleus.P4()->BoostVector();
                auto vec = nu_lep - p4_lep;
                vec.Boost(boost);
                return vec;
              },
              {"gmcrec", "boost"})
      .Define("q0",
              [](const TLorentzVector &q) {
                // return event.event.event->Summary()->KinePtr()->Q2();
                return q.T();
              },
              {"q"})
      .Define("Q2",
              [](const TLorentzVector &q) {
                // return event.event.event->Summary()->KinePtr()->Q2();
                return -q.M2();
              },
              {"q"})
      .Define("visE",
              [](const NtpMCEventRecord &event) {
                double visE{};
                auto nentries = event.event->GetEntriesFast();
                for (int i = 0; i < nentries; i++) {
                  auto part = event.event->Particle(i);
                  if (part->Status() == genie::kIStStableFinalState &&
                      part->Charge() != 0) {
                    // remove mass part for p/n
                    visE += part->E();
                    if (pdg::IsNeutron(part->Pdg()) ||
                        pdg::IsProton(part->Pdg())) {
                      visE -= part->Mass();
                    }
                  }
                }
                return visE;
              },
              {"gmcrec"})
      .Define("nch",
              [](const NtpMCEventRecord &event) {
                int nch{};
                auto nentries = event.event->GetEntriesFast();
                for (int i = 0; i < nentries; i++) {
                  auto part = event.event->Particle(i);
                  if (part->Status() == genie::kIStStableFinalState &&
                      part->Charge() != 0) {
                    nch++;
                  }
                }
                return (double)nch;
              },
              {"gmcrec"})
      // todo: add angle / mag of p_had_system p_max_hadron
      .Define("p_had_system",
              [](const NtpMCEventRecord &event) {
                auto np = event.event->GetEntriesFast();
                TLorentzVector p{};
                for (int i{}; i < np; i++) {
                  auto particle = event.event->Particle(i);
                  if (particle->Status() == kIStStableFinalState &&
                      pdg::IsHadron(particle->Pdg())) {
                    p += *particle->P4();
                  }
                }
                return p;
              },
              {"gmcrec"})
      .Define("angle_p_had_system",
              [](const TLorentzVector &v) { return v.Theta(); },
              {"p_had_system"})
      .Define("mag_p_had_system", [](const TLorentzVector &v) { return v.P(); },
              {"p_had_system"})
      .Define("p_max_hadron",
              [](const NtpMCEventRecord &event) {
                auto np = event.event->GetEntriesFast();
                TLorentzVector p{};
                for (int i{}; i < np; i++) {
                  auto particle = event.event->Particle(i);
                  if (particle->Status() == kIStStableFinalState &&
                      pdg::IsHadron(particle->Pdg())) {
                    if (particle->P4()->P() > p.P()) {
                      p = *particle->P4();
                    }
                  }
                }
                return p;
              },
              {"gmcrec"})
      .Define("angle_p_max_hadron",
              [](const TLorentzVector &v) { return v.Theta(); },
              {"p_max_hadron"})
      .Define("mag_p_max_hadron", [](const TLorentzVector &v) { return v.P(); },
              {"p_max_hadron"})
      .Define("mag_pl_max_hadron",
              [](const TLorentzVector &v) { return v.Pz(); }, {"p_max_hadron"})
      .Define("pt",
              [](const NtpMCEventRecord &event) {
                auto np = event.event->GetEntriesFast();
                TLorentzVector p{};
                for (int i{}; i < np; i++) {
                  auto particle = event.event->Particle(i);
                  if (particle->Status() == kIStStableFinalState &&
                      pdg::IsHadron(particle->Pdg())) {
                    p += *particle->P4();
                  }
                }
                auto pt = p.Pt();
                return pt;
              },
              {"gmcrec"})
      .Define("pl",
              [](const NtpMCEventRecord &event) {
                auto np = event.event->GetEntriesFast();
                TLorentzVector p{};
                for (int i{}; i < np; i++) {
                  auto particle = event.event->Particle(i);
                  if (particle->Status() == kIStStableFinalState &&
                      pdg::IsHadron(particle->Pdg())) {
                    p += *particle->P4();
                  }
                }
                auto pl = p.Pz();
                return pl;
              },
              {"gmcrec"})
      .Define("W",
              [](const NtpMCEventRecord &event) {
                return event.event->Summary()->Kine().W(true);
              },
              {"gmcrec"});
  ;
  // for now we filter out event outside of the binning range
}

// we don't have a proper configuration mechanism yet, so we need to
// hardcode sth.
int main(int argc, char **argv) {
  // ROOT::EnableImplicitMT(4);

  TH1::AddDirectory(false);
  std::string from_id = argc > 1 ? argv[1] : "0000";
  std::string to_id = argc > 2 ? argv[2] : "0001";
  std::string basedir =
      argc > 3 ? argv[3] : "/var/home/yan/neutrino/had_scan_test/had_scan";
  if (basedir.back() != '/')
    basedir += "/";

  auto gevgen_filename = "numu_on_p.ghep.root";
  // auto gevgen_filename = "gntp.70000.ghep.root";

  std::string input_from_file =
      basedir + "scan/" + from_id +
      "/3.04.02-routine_validation_01-hadronization/" + gevgen_filename;
  std::string input_to_file = basedir + "scan/" + to_id +
                              "/3.04.02-routine_validation_01-hadronization/" +
                              gevgen_filename;
  auto spline_path_from =
      basedir + "scan/" + from_id +
      "/3.04.02-routine_validation_01-xsec_vA/total_xsec.root";
  auto spline_path_to =
      basedir + "scan/" + to_id +
      "/3.04.02-routine_validation_01-xsec_vA/total_xsec.root";
  auto ipol_path = basedir + "ipol_test.dat";
  auto binningxml = basedir + "binning.xml";

  auto input_from_tree =
      ROOT::RDataFrame{"gtree", input_from_file}; // read the tree from the file
  ROOT::RDF::Experimental::AddProgressBar(input_from_tree);
  auto observable_splines = std::make_unique<GReWeightProfessor>("");
  observable_splines->ReadComparisonXML(binningxml, ipol_path);
  LOG("Test", pNOTICE) << "Read Professor Reweight success";

  auto from = fromdat(basedir + "scan/" + from_id + "/params.dat");
  auto to = fromdat(basedir + "scan/" + to_id + "/params.dat");

  observable_splines->SetSystematic(to, from);
  LOG("Test", pNOTICE) << "Set Systematic success";

  LOG("Test", pNOTICE) << "Testing RDataFrame::Define to calculate weight\n";
  auto tree_to = ROOT::RDataFrame{"gtree", input_to_file};
  ROOT::RDF::Experimental::AddProgressBar(tree_to);

  auto ds_ref = common_def(tree_to);

  auto tree_rew =
      common_def(input_from_tree)
          .Define(
              "weight",
              [&](NtpMCEventRecord &event) {
                auto weight = observable_splines->CalcWeight(*(event.event));
                // if (abs(weight - 1.) > 1e-3)
                //   std::cout << "At lease we have " << weight
                //   << std::endl;
                // event.Clear();
                return (weight < 0 || std::isnan(weight) || std::isinf(weight))
                           ? 0.
                           : weight;
                // return weight;
              },
              {"gmcrec"});

  auto modelmap = std::to_array<std::tuple<std::string, ROOT::RDF::TH1DModel>>(
      {{"muonp", {"muonp", ";p_{#mu};", 10, 0, 10.}},
       {"muon_theta", {"muon_theta", ";#theta_{#mu};", 20, 0., 2.}},
       {"Q2", {"Q2", ";Q^{2};", 30, 0., 6.}},
       {"q0", {"q0", ";q^{0};", 30, 0., 6.}},
       {"visE", {"visE", ";E_{vis};", 10, 0., 100.}},
       {"visE", {"visE_low", ";E_{vis};", 10, 0., 10.}},
       {"nch", {"nch", ";nch;", 19, -0.5, 18.5}},
       {"pt", {"pt", ";pt;", 11, 0., 10.}},
       {"pl", {"pl", ";pl;", 11, 0., 10.}},
       {"W", {"W", ";W;", 20, 0., 30.}},
       {"W", {"W_low", ";W;", 10, 0., 6.}},
       {"W", {"W_3low", ";W;", 10, 0., 3.}},
       {"angle_p_max_hadron",
        {"angle_p_max_hadron", ";#theta_{p_{max}};", 10, 0, 1.5}},
       {"mag_p_max_hadron", {"mag_p_max_hadron", ";|p_{max}|;", 10, 0., 10.}},
       {"angle_p_had_system",
        {"angle_p_had_system", ";#theta_{p_{had}};", 10, 0, 1.5}},
       {"mag_p_had_system", {"mag_p_had_system", ";|p_{had}|;", 10, 0., 10.}},
       {"mag_pl_max_hadron",
        {"mag_pl_max_hadron", ";|p_{had}|;", 10, 0., 10.}}});

  // auto

  auto merge_cut_lists = [merge_cut_func = [](auto &&...args)
                              -> std::function<bool(NtpMCEventRecord &)> {
    return [args...](NtpMCEventRecord &event) { return (args(event) && ...); };
  }](auto &&lista, auto &&listb) {
    return std::views::transform(
               std::views::cartesian_product(lista, listb),
               [merge_cut_func](auto &&cuts) {
                 auto [cuta, cutb] = cuts;
                 auto [namea, cutfunc_a] = cuta;
                 auto [nameb, cutfunc_b] = cutb;
                 auto newcut = merge_cut_func(cutfunc_a, cutfunc_b);
                 return std::make_tuple(namea + nameb, newcut);
               }) |
           std::ranges::to<std::vector>();
  };

  auto W_cuts = std::to_array<
      std::tuple<std::string, std::function<bool(NtpMCEventRecord &)>>>(
      {{"LowW",
        [](const NtpMCEventRecord &event) -> bool {
          double W = event.event->Summary()->Kine().W(true);
          return W < 2.3;
        }},
       {"MidW",
        [](const NtpMCEventRecord &event) -> bool {
          double W = event.event->Summary()->Kine().W(true);
          return W > 2.3 && W < 3;
        }},
       {"HiW", [](const NtpMCEventRecord &event) -> bool {
          double W = event.event->Summary()->Kine().W(true);
          return W > 3;
        }}});

  auto to_npi = [](const NtpMCEventRecord &event) {
    size_t npi{};
    auto np = event.event->GetEntriesFast();
    for (int i{}; i < np; i++) {
      auto particle = event.event->Particle(i);
      if (particle->Status() == kIStStableFinalState &&
          pdg::IsPion(particle->Pdg())) {
        npi++;
      }
    }
    return npi;
  };

  auto multi_cuts = std::to_array<
      std::tuple<std::string, std::function<bool(NtpMCEventRecord &)>>>(
      {{"0pi",
        [&](const NtpMCEventRecord &event) -> bool {
          return to_npi(event) == 0;
        }},
       {"1pi",
        [&](const NtpMCEventRecord &event) -> bool {
          return to_npi(event) == 1;
        }},
       {"2pi",
        [&](const NtpMCEventRecord &event) -> bool {
          return to_npi(event) == 2;
        }},
       {"Mpi", [&](const NtpMCEventRecord &event) -> bool {
          return to_npi(event) >= 3;
        }}});

  auto channelcuts = merge_cut_lists(W_cuts, multi_cuts);
  channelcuts.insert(channelcuts.end(), W_cuts.begin(), W_cuts.end());

  //   ,
  //  {"LowW_0pi",
  //   [](const NtpMCEventRecord &event) -> bool {
  //     double W = event.event->Summary()->Kine().W(true);
  //     size_t npi{};
  //     auto np = event.event->GetEntriesFast();
  //     for (int i{}; i < np; i++) {
  //       auto particle = event.event->Particle(i);
  //       if (particle->Status() == kIStStableFinalState &&
  //           pdg::IsPion(particle->Pdg())) {
  //         npi++;
  //       }
  //     }
  //     return W < 3 && npi == 0;
  //   }},
  //  {"LowW_1pi",
  //   [](const NtpMCEventRecord &event) -> bool {
  //     double W = event.event->Summary()->Kine().W(true);
  //     size_t npi{};
  //     auto np = event.event->GetEntriesFast();
  //     for (int i{}; i < np; i++) {
  //       auto particle = event.event->Particle(i);
  //       if (particle->Status() == kIStStableFinalState &&
  //           pdg::IsPion(particle->Pdg())) {
  //         npi++;
  //       }
  //     }
  //     return W < 3 && npi == 1;
  //   }},
  //  {"LowW_2pi",
  //   [](const NtpMCEventRecord &event) -> bool {
  //     double W = event.event->Summary()->Kine().W(true);
  //     size_t npi{};
  //     auto np = event.event->GetEntriesFast();
  //     for (int i{}; i < np; i++) {
  //       auto particle = event.event->Particle(i);
  //       if (particle->Status() == kIStStableFinalState &&
  //           pdg::IsPion(particle->Pdg())) {
  //         npi++;
  //       }
  //     }
  //     return W < 3 && npi == 2;
  //   }},
  //  {"LowW_Mpi", [](const NtpMCEventRecord &event) -> bool {
  //     double W = event.event->Summary()->Kine().W(true);
  //     size_t npi{};
  //     auto np = event.event->GetEntriesFast();
  //     for (int i{}; i < np; i++) {
  //       auto particle = event.event->Particle(i);
  //       if (particle->Status() == kIStStableFinalState &&
  //           pdg::IsPion(particle->Pdg())) {
  //         npi++;
  //       }
  //     }
  //     return W < 3 && npi >= 3;
  //   }}}

  auto hists_ref =
      modelmap |
      std::views::transform(
          [&](const std::tuple<std::string, ROOT::RDF::TH1DModel> &obj) {
            auto [var, model] = obj;
            auto plot_cut =
                channelcuts |
                std::views::transform(
                    [&](const std::tuple<
                        std::string, std::function<bool(NtpMCEventRecord &)>>
                            &obj) {
                      auto &[cutname, cut] = obj;
                      auto thismodel = model;
                      thismodel.fName += cutname;
                      return ds_ref.Filter(cut, {"gmcrec"})
                          .Histo1D<double>(thismodel, var);
                    }) |
                std::ranges::to<std::vector>();
            return std::make_tuple(ds_ref.Histo1D<double>(model, var),
                                   plot_cut);
          }) |
      std::ranges::to<std::vector>();

  auto hists_rew =
      modelmap |
      std::views::transform(
          [&](const std::tuple<std::string, ROOT::RDF::TH1DModel> &obj) {
            auto [var, model] = obj;
            auto plot_cut =
                channelcuts |
                std::views::transform(
                    [&](const std::tuple<
                        std::string, std::function<bool(NtpMCEventRecord &)>>
                            &obj) {
                      auto &[cutname, cut] = obj;
                      auto thismodel = model;
                      thismodel.fName += cutname;
                      return tree_rew.Filter(cut, {"gmcrec"})
                          .Histo1D<double, double>(thismodel, var, "weight");
                    }) |
                std::ranges::to<std::vector>();
            return std::make_tuple(
                tree_rew.Histo1D<double, double>(model, var, "weight"),
                plot_cut);
          }) |
      std::ranges::to<std::vector>();

  auto hists_unrew =
      modelmap |
      std::views::transform(
          [&](const std::tuple<std::string, ROOT::RDF::TH1DModel> &obj) {
            auto [var, model] = obj;
            auto plot_cut =
                channelcuts |
                std::views::transform(
                    [&](const std::tuple<
                        std::string, std::function<bool(NtpMCEventRecord &)>>
                            &obj) {
                      auto &[cutname, cut] = obj;
                      auto thismodel = model;
                      thismodel.fName += cutname;
                      return tree_rew.Filter(cut, {"gmcrec"})
                          .Histo1D<double>(thismodel, var);
                    }) |
                std::ranges::to<std::vector>();
            return std::make_tuple(tree_rew.Histo1D<double>(model, var),
                                   plot_cut);
          }) |
      std::ranges::to<std::vector>();

  std::cout << "Running RDF\n\n" << std::flush;
  ROOT::RDF::RunGraphs({
      std::get<0>(hists_ref[0]),
      std::get<0>(hists_rew[0]),
      std::get<0>(hists_unrew[0]),
  });

  auto do_plot = [&](TH1D *hist_rew, TH1D *hist_reference, TH1D *hist_unrew) {
    std::string ytitle = "Count (Count * Weight)";

    auto chi2 = hist_reference->Chi2Test(hist_rew, "CHI2/NDF");

    auto c1 = std::make_unique<TCanvas>();
    c1->Draw();

    constexpr double fair_share = 0.7;
    auto pad1 = std::make_unique<TPad>("pad1", "", 0, 1 - fair_share, 1, 1);
    auto pad2 = std::make_unique<TPad>("pad2", "", 0, 0, 1, 1 - fair_share);
    pad1->SetBottomMargin(0);
    pad1->SetTopMargin(pad1->GetTopMargin() / 2.);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(pad2->GetBottomMargin() * 2.);

    // c1->cd(1);
    pad1->cd();
    auto max = std::max({hist_reference->GetMaximum(), hist_rew->GetMaximum(),
                         hist_unrew->GetMaximum()});
    gStyle->SetOptStat(0);
    auto legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
    double ytitleoffset = 1.2;
    hist_rew->SetLineColor(kBlue);
    hist_rew->SetMaximum(max * 1.1);
    hist_rew->GetYaxis()->SetTitle(ytitle.c_str());
    hist_rew->GetYaxis()->SetTitleOffset(ytitleoffset);
    hist_rew->Draw("hist E");
    legend->AddEntry(hist_rew, "reweighted", "l");
    hist_unrew->SetLineColor(kGreen);
    hist_unrew->SetMaximum(max * 1.1);
    hist_unrew->SetLineStyle(kDashed);
    // hist_unrew->SetLineStyle(kDotted);
    hist_unrew->GetYaxis()->SetTitle(ytitle.c_str());
    hist_unrew->GetYaxis()->SetTitleOffset(ytitleoffset);
    hist_unrew->Draw("hist E same");
    legend->AddEntry(hist_unrew, "unreweighted", "l");
    hist_reference->SetMaximum(max * 1.1);
    hist_reference->SetLineColor(kRed);
    // hist_reference->SetLineStyle(kDotted);
    hist_reference->GetYaxis()->SetTitle(ytitle.c_str());
    hist_reference->GetYaxis()->SetTitleOffset(ytitleoffset);
    hist_reference->SetLineStyle(kDashed);
    hist_reference->Draw("hist E SAME");
    legend->AddEntry(hist_reference, "reference", "l");
    legend->SetHeader(("chi2/ndf: " + std::to_string(chi2)).c_str());
    legend->Draw("SAME");
    // c1->cd(2);
    pad2->cd();
    auto hist_diff_rew_ref =
        std::unique_ptr<TH1D>{static_cast<TH1D *>(hist_rew->Clone())};
    hist_diff_rew_ref->Add(hist_reference, -1);
    // hist_diff_rew_ref: Rew - Ref

    // auto hist_diff_unrew_ref =
    //     std::unique_ptr<TH1D>{static_cast<TH1D *>(hist_unrew->Clone())};
    // hist_diff_unrew_ref->Add(hist_reference, -1);
    // (hist_rew - hist_reference) / (hist_unrew - hist_reference);
    hist_diff_rew_ref->Divide(hist_reference);
    // hist_diff_rew_ref->SetMaximum();
    // hist_diff_rew_ref->Set
    // for (int bin_id = 1; bin_id < hist_diff_rew_ref->GetNbinsX(); bin_id++) {
    //   double sigma_rew = hist_rew->GetBinError(bin_id);
    //   double sigma_unrew = hist_unrew->GetBinError(bin_id);
    //   double sigma_ref = hist_reference->GetBinError(bin_id);
    //   double rew = hist_rew->GetBinContent(bin_id);
    //   double unrew = hist_unrew->GetBinContent(bin_id);
    //   double ref = hist_reference->GetBinContent(bin_id);
    //   using std::sqrt, std::pow;
    //   double new_err = sqrt(
    //       (pow(1 / (ref - unrew), 2) * pow(sigma_rew, 2)) +
    //       (pow((rew - ref) / pow(ref - unrew, 2), 2) * pow(sigma_unrew, 2)) +
    //       (pow(-(rew - unrew) / pow(ref - unrew, 2), 2) * pow(sigma_ref,
    //       2)));
    //   if (rew == 0 || unrew == 0 || ref == 0) {
    //     new_err = 0;
    //   }
    //   hist_diff_rew_ref->SetBinError(bin_id, new_err);
    // }
    auto text_size = hist_diff_rew_ref->GetYaxis()->GetLabelSize();
    auto title_offset = hist_reference->GetYaxis()->GetTitleOffset();

    auto text_size_scaled = text_size * fair_share / (1 - fair_share);
    auto title_offset_scaled = title_offset / (fair_share / (1 - fair_share));
    hist_diff_rew_ref->GetXaxis()->SetLabelSize(text_size_scaled);
    hist_diff_rew_ref->GetXaxis()->SetTitleSize(text_size_scaled);
    hist_diff_rew_ref->GetYaxis()->SetLabelSize(text_size_scaled);
    hist_diff_rew_ref->GetYaxis()->SetTitleSize(text_size_scaled);

    hist_diff_rew_ref->GetYaxis()->SetTitleOffset(title_offset_scaled);
    hist_diff_rew_ref->GetYaxis()->SetTitle("rew - ref / ref");
    hist_diff_rew_ref->GetYaxis()->CenterTitle(true);

    // hist_diff_rew_ref->SetMaximum(hist_diff_rew_ref->GetMaximum() + 0.5);
    // hist_diff_rew_ref->SetMinimum(hist_diff_rew_ref->GetMinimum() - 0.5);

    hist_diff_rew_ref->Draw("histE1");
    c1->cd();
    pad1->Draw();
    pad2->Draw();
    c1->Draw();
    // std::string filename =
    //     varname + "_" + (do_normalize ? "normalized" : "unnormalized") + cut;
    std::string filename = hist_rew->GetName();
    c1->SaveAs((filename + ".pdf").c_str());
    c1->SaveAs((filename + ".eps").c_str());
  };

  for (auto &&[var, ent_ref, ent_rew, ent_unrew] : std::views::zip(
           modelmap | std::views::keys, hists_ref, hists_rew, hists_unrew)) {
    auto [hist_ref, cut_plot_ref] = ent_ref;
    auto [hist_rew, cut_plot_rew] = ent_rew;
    auto [hist_unrew, cut_plot_unrew] = ent_unrew;
    std::cout << "Plotting " << var << '\n';
    do_plot(hist_rew.GetPtr(), hist_ref.GetPtr(), hist_unrew.GetPtr());
    for (auto &&[cutname, hist_ref, hist_rew, hist_unrew] :
         std::views::zip(channelcuts | std::views::keys, cut_plot_ref,
                         cut_plot_rew, cut_plot_unrew)) {
      std::cout << "Plotting " << var << " with cut: " << cutname << '\n';
      do_plot(hist_rew.GetPtr(), hist_ref.GetPtr(), hist_unrew.GetPtr());
    }
  }

  // tree_rew.Snapshot("out", "out.root", {"weight", "Q2"});

  // // do_plot("muonp", true);
  // // do_plot("muonp", false);
  // // do_plot("muon_theta", true);
  // // do_plot("muon_theta", false);
  // // do_plot("Q2", true);
  // // do_plot("Q2", false);
  // for (auto &&var : varnames) {
  //   // do_plot(var, false);
  //   do_plot(var, false);
  //   do_plot(var, false, "LowW");
  //   do_plot(var, false, "HiW");
  // }

  return 0;
}