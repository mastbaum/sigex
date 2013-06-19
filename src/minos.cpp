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
#include <TGraph.h>
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
Fit* run_experiment(std::vector<Signal> signals, std::string signal_name, float live_time,
                    Range<float> e_range, Range<float> r_range,
                    double& norm, double& err) {
  // make a dataset with zero signal
  FakeDataGenerator gen(signals, live_time, e_range, r_range);
  std::map<std::string, double> mnorm = {{ signal_name, 0 }};
  TNtuple* nt = gen.make_dataset(false, &mnorm);
  Dataset data("data", nt);

  // best fit
  Fit* fit = new Fit(signals, data);
  TMinuit* minuit = (*fit)(e_range, r_range, live_time, nullptr, nullptr, true);  // compute minos errors

  double lfit, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat(lfit, edm, errdef, nvpar, nparx, icstat);
  std::cout << "Best fit:" << std::endl;
  minuit->mnprin(4, lfit);

  int eflag;
  minuit->mnexcm("SHOW COR", static_cast<double*>(nullptr), 0, eflag);

  double eplus, eminus, eparabolic, globcc;
  for (size_t i=0; i<signals.size(); i++) {
    if (signals[i].name == signal_name) {
      double signal_par = 0;
      minuit->GetParameter(i, signal_par, eparabolic);
      norm = signal_par * signals[i].nexpected;
      minuit->mnerrs(i, eplus, eminus, eparabolic, globcc);
      err = eplus * signals[i].nexpected;
      break;
    }
  }

  //double* args = new double[4];
  //args[0] = 1;
  //args[1] = 9;
  //args[2] = 1.645;
  //args[3] = 50;
  //minuit->mnexcm("CONTOUR", args, 4, eflag);

  if (std::string(minuit->fCstatu) == "SUCCESSFUL") {
    return fit;
  }
  else {
    std::cerr << "Warning: MINUIT failed with status " << minuit->fCstatu << std::endl;
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
  std::map<double, double> quantiles = {{    0,     0 },
                                        { 6800, 1.000 },
                                        { 9000, 1.645 },
                                        { 9500, 1.960 },
                                        { 9900, 2.576 },
                                        { 9990, 3.291 },
                                        { 9999, 3.891 }};
  int confidence = static_cast<int>(fc.confidence * 10000);
  double quantile = quantiles.at(confidence);

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
    Fit* fit = run_experiment(fc.signals, fc.signal_name, fc.live_time, fc.e_range, fc.r_range, norm, err);
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


  // plots
  gStyle->SetOptStat(0);
  std::vector<int> colors = {2, 3, 4, 6, 7, 8, 9, 11, 29, 5, 30, 34, 38, 40, 42, 45, 49};

  TCanvas ce;
  ce.SetLogy();
  TCanvas cr;
  cr.SetLogy();
  TLegend stan(0.85, 0.35, 0.98, 0.88);
  stan.SetBorderSize(1);
  stan.SetFillColor(0);

  TH1F* hesum = nullptr;
  TH1F* hesum_ns = nullptr;
  TH1F* hrsum = nullptr;
  TH1F* hrsum_ns = nullptr;

  for (size_t i=0; i<fc.signals.size(); i++) {
    TH1F* he = nullptr;
    TH1F* hr = nullptr;
    TH1* ht = fc.signals[i].histogram;
    if (ht->IsA() == TH2F::Class()) {
      TH2F* h2 = (TH2F*)(ht);
      int first_bin = h2->GetXaxis()->FindBin(fc.r_range.min);
      int last_bin = h2->GetXaxis()->FindBin(fc.r_range.max);
      TH1D* tt = h2->ProjectionY("pdf_proj_y", first_bin, last_bin);
      TH1F temp;
      tt->Copy(temp);  // convert TH1D to TH1F
      he = (TH1F*) temp.Clone("he");

      first_bin = ht->GetYaxis()->FindBin(fc.e_range.min);
      last_bin = ht->GetYaxis()->FindBin(fc.e_range.max);
      TH1D* temp2 = dynamic_cast<TH2F*>(ht)->ProjectionX("hr", first_bin, last_bin);
      TH1F temp3;
      temp2->Copy(temp3);
      hr = new TH1F(temp3);
    }
    else {
      he = (TH1F*)(ht);
    }

    double p = (fc.signals[i].name == fc.signal_name ? meds : fc.signals[i].nexpected);
//    std::cout << fc.signals[i].name << ": " << p << " " << integral << " " << he->Integral() << std::endl;
    std::cout << ">>> " << fc.signals[i].name << " " << p << " " << he->Integral() << std::endl;
    he->Scale(p/he->Integral());
    if (hr) {
      hr->Scale(p/hr->Integral());
    }

    if (hesum == nullptr) {
      if (hr) {
        hrsum = (TH1F*)hr->Clone("hrsum");
        hrsum->Reset();
        hrsum->SetLineColor(kBlack);
        hrsum->SetLineWidth(2);

        hrsum_ns = (TH1F*)hr->Clone("hrsum_ns");
        hrsum_ns->Reset();
        hrsum_ns->SetLineColor(kBlack);
        hrsum_ns->SetLineWidth(2);
        hrsum_ns->SetLineStyle(2);
      }

      hesum = (TH1F*)he->Clone("hesum");
      hesum->Reset();
      hesum->SetLineColor(kBlack);
      hesum->SetLineWidth(2);
      stan.AddEntry(hesum, "Sum");

      hesum_ns = (TH1F*)he->Clone("hesum_ns");
      hesum_ns->Reset();
      hesum_ns->SetLineColor(kBlack);
      hesum_ns->SetLineWidth(2);
      hesum_ns->SetLineStyle(2);
      stan.AddEntry(hesum_ns, "Sum, no signal");
    }

    if (hr) {
      hrsum->Add(hr);
      if (fc.signals[i].name != fc.signal_name) {
        hrsum_ns->Add(hr);
      }
    }

    hesum->Add(he);
    if (fc.signals[i].name != fc.signal_name) {
      hesum_ns->Add(he);
    }

    if (hr) {
      hr->GetXaxis()->SetRangeUser(fc.r_range.min, fc.r_range.max);
      hr->SetAxisRange(5e-2, 1e3, "Y"); //1e-4, 5, "Y");
      hr->SetLineColor(colors[i%colors.size()]);
      hr->SetLineWidth(2);
    }

    he->GetXaxis()->SetRangeUser(2, 3.2);
    he->SetAxisRange(1e-2, 5e2, "Y");
    char cl[100];
    snprintf(cl, 100, "%1.2f%% CL #hat{S} = %1.2f (%i experiments, %1.0f y)",
             fc.confidence*100, meds, fc.mc_trials, fc.live_time);
    he->SetTitle(cl);
    he->SetLineColor(colors[i%colors.size()]);
    he->SetLineWidth(2);
    stan.AddEntry(he, fc.signals[i].title.c_str());

    ce.cd();
    he->Draw(i==0 ? "hist" : "hist same");

    if (hr) {
      cr.cd();
      hr->Draw(i==0 ? "hist" : "hist same");
    }
  }

  ce.cd();
  hesum->Draw("hist same");
  hesum_ns->Draw("hist same");

  FakeDataGenerator gen(fc.signals, fc.live_time, fc.e_range, fc.r_range);
  std::map<std::string, double> mnorm = {{ fc.signal_name, 0 }};
  TNtuple* nt = gen.make_dataset(false, &mnorm);
  TH1F* hfd = (TH1F*) hesum->Clone("hfd");
  hfd->Reset();
  hfd->SetMarkerStyle(20);
  nt->Draw("e>>hfd","","same");

  stan.Draw();
  ce.SaveAs((fc.output_file + "_e.pdf").c_str());

  if (hrsum) {
    cr.cd();
    hrsum->Draw("hist same");
    hrsum_ns->Draw("hist same");
    stan.Draw();
    cr.SaveAs((fc.output_file + "_r.pdf").c_str());
  }

/*
  // energy
  TCanvas c2;
  c2.SetLogy();

  TH1F* hsum = h.;
  TH1F* hsum_ns = nullptr;
  TLegend stan(0.75, 0.35, 0.88, 0.88);
  stan.SetBorderSize(1);
  stan.SetFillColor(0);
  bool first = true;
  size_t icolor = 0;
  for (size_t i=0; i<fc.signals.size(); i++) {
    TH1* ht = fc.signals[i].histogram;
    TH1F h;
    if (ht->IsA() == TH2F::Class()) {
      h.Copy(*FitConfig::project1d(dynamic_cast<TH2F*>(ht), &fc.r_range));
    }
    else if (h.IsA() == TH1F::Class()) {
      h = *((TH1F*)(ht));
    }

    double p = (fc.signals[i].name == fc.signal_name ? meds : fc.signals[i].rate);
    int x1 = h.FindBin(fc.e_range.min);
    int x2 = h.FindBin(fc.e_range.max);
    double ii = h.Integral(x1, x2);
    h.Scale(p/ii);

    h.SetAxisRange(1e-1, 1e4, "Y");
    char binsize[20];
    snprintf(binsize, 20, "%1.0f", h.GetXaxis()->GetBinWidth(1) * 1000);
    h.SetYTitle(("Counts per " + std::string(binsize) + " keV bin").c_str());
    h.GetXaxis()->SetRangeUser(2, 3.2);
    char cl[100];
    snprintf(cl, 100, "%1.2f%% CL #hat{S} = %1.2f (%i experiments, %1.0f y)",
             fc.confidence*100, meds, fc.mc_trials, fc.live_time);
    h.SetTitle(cl);
    h.SetXTitle("Energy (MeV)");

    if (hsum == nullptr) {
      hsum = (TH1F*) h.Clone("hsum");
      hsum->Reset();
      hsum->SetLineColor(1);
      hsum->SetLineWidth(2);
      stan.AddEntry(hsum, "Sum");
    }

    if (hsum_ns == nullptr) {
      hsum_ns = (TH1F*) h.Clone("hsum_ns");
      hsum_ns->Reset();
      hsum_ns->SetLineColor(1);
      hsum_ns->SetLineWidth(2);
      hsum_ns->SetLineStyle(2);
      stan.AddEntry(hsum_ns, "Sum, no signal");
    }

    h.SetLineColor(colors[icolor%colors.size()]);
    h.SetLineWidth(2);
    h.Draw(first ? "": "same");
    stan.AddEntry(&h, fc.signals[i].title.c_str());

    hsum->Add(&h);
    if (fc.signals[i].name != fc.signal_name) {
      hsum_ns->Add(&h);
    }

    icolor++;
    first = false;
  }

  hsum->Draw("same");
  hsum_ns->Draw("same");
  stan.Draw();
  c2.Update();
  c2.SaveAs((fc.output_file + ".C").c_str());
  c2.SaveAs((fc.output_file + ".pdf").c_str());
*/
  return 0;
}

