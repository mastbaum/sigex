// estimate sensitivity with a single ML fit, via MINOS errors

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <assert.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TMatrixDSym.h>
#include <TF1.h>

#include "util.h"
#include "fit.h"

#ifndef VERBOSE
#define VERBOSE false
#endif

// run an "experiment:" generate a dataset sampled from histograms and do an
// unbinned ML fit. returns the fit object if successful.
Fit* run_experiment(std::vector<Signal> signals, std::string signal_name,
                    Range<float> e_range, Range<float> r_range,
                    double& norm, double& err) {
  // make a dataset with zero signal
  FakeDataGenerator gen(signals, e_range, r_range);
  std::map<std::string, double> mnorm = {{ signal_name, 0 }};
  TNtuple* nt = gen.make_dataset(true, &mnorm);
  Dataset data("data", nt);

  // best fit
  Fit* fit = new Fit(signals, data);
  const TMinuit* minuit = (*fit)(e_range, nullptr, nullptr, true);  // compute minos errors

  double lfit, edm, errdef;
  int nvpar, nparx, icstat;
  const_cast<TMinuit*>(minuit)->mnstat(lfit, edm, errdef, nvpar, nparx, icstat);
  std::cout << "Best fit:" << std::endl;
  const_cast<TMinuit*>(minuit)->mnprin(4, lfit);

  int eflag;
  const_cast<TMinuit*>(minuit)->mnexcm("SHOW COR", static_cast<double*>(nullptr), 0, eflag);

  double eplus, eminus, eparabolic, globcc;
  for (size_t i=0; i<signals.size(); i++) {
    if (signals[i].name == signal_name) {
      double signal_par = 0;
      minuit->GetParameter(i, signal_par, eparabolic);
      norm = signal_par * signals[i].rate;
      const_cast<TMinuit*>(minuit)->mnerrs(i, eplus, eminus, eparabolic, globcc);
      err = eplus * signals[i].rate;
      break;
    }
  }

  if (std::string(minuit->fCstatu) == "SUCCESSFUL") {
    return fit;
  }
  else {
    delete fit;
    return static_cast<Fit*>(nullptr);
  }
}


int main(int argc, char* argv[]) {
  gRandom->SetSeed(0);

  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " config.json" << std::endl;
    return 1;
  }

  // load configuration
  const std::string config_filename = std::string(argv[1]);
  FitConfig fc(config_filename);
  fc.print();

  // convert central-value Gaussian CL (x 10000) to quantile
  // cf. Cowan p. 125
  // todo: replace with an error function
  std::map<double, double> quantiles = {{ 6800, 1.000 },
                                        { 9000, 1.645 },
                                        { 9500, 1.960 },
                                        { 9900, 2.576 },
                                        { 9990, 3.291 },
                                        { 9999, 3.891 }};
  int confidence = static_cast<int>(fc.confidence * 10000);
  double quantile = quantiles[confidence];

  // run an ensemble of fake experiments
  TH1F hs("hs", "#hat{S}", 550, -50, 500);
  TH1F hs_shift("hs_shift", "#hat{S}, Shift method", 550, -50, 500);
  std::vector<double> hats;
  std::vector<double> hats_shift;

  TH1F means("means", "means", 10000, -5000, 5000);
  TH1F widths("widths", "widths", 10000, 0, 10000);

  for (size_t i=0; i<fc.mc_trials; i++) {
    std::cout << "=== EXPERIMENT " << i << " ===" << std::endl;
    double norm = 0;
    double err = 0;
    Fit* fit = run_experiment(fc.signals, fc.signal_name, fc.e_range, fc.r_range, norm, err);
    if (fit != static_cast<Fit*>(nullptr)) {
      std::cout << " Best fit: " << norm << ", positive error " << err << std::endl;

      // build up some histograms for the nasim method
      means.Fill(norm);
      widths.Fill(err);

      // just use minos results directly
      // some "lucky" experiments will set negative limits
      double s = norm + quantile * err;
      hs.Fill(s);
      hats.push_back(s);

      // shift negative means to zero since everyone is doing it.
      // cf. Statistical Data Analysis, Glen Cowan, p. 137
      double s_shift = std::max(norm, 0.0) + quantile * err;
      hs_shift.Fill(s_shift);
      hats_shift.push_back(s_shift);
    }
    else {
      std::cout << " Skipping fit due to MINUIT errors" << std::endl;
    }
    delete fit;
  }

  TCanvas c1;
  hs.DrawNormalized();
  hs_shift.SetLineColor(kBlue);
  hs_shift.DrawNormalized("same");
  c1.SaveAs("hats.pdf");
  c1.SaveAs("hats.C");

  double meds = median(hats);
  std::cout << "Average limit: " << meds << std::endl;

  double meds_shift = median(hats_shift);
  std::cout << "Average limit (shift method): " << meds_shift << std::endl;

  means.Fit("gaus");
  widths.Fit("gaus");
  TF1* f = means.GetFunction("gaus");
  double meds_msum = 0;
  if (f != nullptr) {
    meds_msum = f->GetParameter(1) + quantile * widths.GetFunction("gaus")->GetParameter(1);
  }
  std::cout << "Average limit (sum of means method): " << meds_msum << std::endl;


  // plot
  //FakeDataGenerator gen(fc.signals, fc.e_range, fc.r_range);
  //std::map<std::string, double> mnorm = {{ fc.signal_name, 0 }};
  //TNtuple* nt = gen.make_dataset(true, &mnorm);
  gStyle->SetOptStat(0);
  std::vector<int> colors = {2, 3, 4, 6, 7, 8, 9, 11, 29, 5, 30, 34, 38, 40, 42, 45, 49};
  TCanvas c2;
  TH1F* hsum = nullptr;
  TH1F* hsum_ns = nullptr;
  TH1F* hdata = nullptr;
  c2.SetLogy();
  TLegend stan(0.65, 0.55, 0.88, 0.88);
  stan.SetBorderSize(0);
  stan.SetFillColor(0);
  bool first = true;
  size_t icolor = 0;
  for (size_t i=0; i<fc.signals.size(); i++) {
    double p = (fc.signals[i].name == fc.signal_name ? meds : fc.signals[i].rate);
    int x1 = fc.signals[i].histogram->FindBin(fc.e_range.min);
    int x2 = fc.signals[i].histogram->FindBin(fc.e_range.max);
    double ii = fc.signals[i].histogram->Integral(x1, x2);
    fc.signals[i].histogram->Scale(p/ii);

    fc.signals[i].histogram->SetAxisRange(1e-1, 1e3, "Y");
    char binsize[20];
    snprintf(binsize, 20, "%1.0f", fc.signals[i].histogram->GetXaxis()->GetBinWidth(1) * 1000);
    fc.signals[i].histogram->SetYTitle(("Counts per " + std::string(binsize) + " keV bin").c_str());
    fc.signals[i].histogram->GetXaxis()->SetRangeUser(2, 3.2);
    char cl[100];
    snprintf(cl, 100, "%1.2f%% CL #hat{S} = %1.2f (%i experiments, %1.0f y)",
             fc.confidence*100, meds, fc.mc_trials, fc.live_time);
    fc.signals[i].histogram->SetTitle(cl);
    fc.signals[i].histogram->SetXTitle("Energy (MeV)");

    if (hdata == nullptr) {
      hdata = (TH1F*) fc.signals[i].histogram->Clone("hdata");
      hdata->SetLineColor(4);
      hdata->SetLineWidth(2);
      //stan.AddEntry(hdata, "Sample MC");
      //nt->Draw("e>>hdata");
      //hdata->Draw(first ? "" : "same");
      //first = false;
    }

    if (hsum == nullptr) {
      hsum = (TH1F*) fc.signals[i].histogram->Clone("hsum");
      hsum->SetLineColor(1);
      hsum->SetLineWidth(2);
      stan.AddEntry(hsum, "Sum");
    }

    if (hsum_ns == nullptr) {
      hsum_ns = (TH1F*) fc.signals[i].histogram->Clone("hsum_ns");
      hsum_ns->SetLineColor(1);
      hsum_ns->SetLineWidth(2);
      hsum_ns->SetLineStyle(2);
      stan.AddEntry(hsum_ns, "Sum, no signal");
    }

    fc.signals[i].histogram->SetLineColor(colors[icolor]);
    fc.signals[i].histogram->SetLineWidth(2);
    fc.signals[i].histogram->Draw(first ? "": "same");
    stan.AddEntry(fc.signals[i].histogram, fc.signals[i].title.c_str());

    hsum->Add(fc.signals[i].histogram);
    if (fc.signals[i].name != fc.signal_name) {
      hsum_ns->Add(fc.signals[i].histogram);
    }

    icolor++;
    first = false;
  }

  hsum->Draw("same");
  hsum_ns->Draw("same");
  stan.Draw();
  c2.Update();
  c2.SaveAs("fit_spectrum_b.C");
  c2.SaveAs("fit_spectrum_b.pdf");

  return 0;
}

