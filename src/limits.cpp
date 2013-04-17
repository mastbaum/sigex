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
#include <TLegend.h>
#include <TRandom.h>
#include <TMinuit.h>

#include "util.h"
#include "fit.h"

#ifndef VERBOSE
#define VERBOSE false
#endif

// run an experiment, get an upper limit
float get_ul(std::vector<Signal> _signals, std::string signal_name,
             float confidence, FakeDataGenerator& gen, unsigned nmc, TH2F* hr,
             Range<float> e_range, Range<float> r_range) {
  std::vector<Signal> signals = _signals;

  // make experiment data
  std::map<std::string, double> snorm = {{signal_name, 0}};
  TNtuple* nt = gen.make_dataset(true, &snorm);
  Dataset data("data", nt);

  // best fit
  Fit* fit = new Fit(signals, data);
  TMinuit* minuit = (*fit)(e_range, r_range);
  double lfit, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat(lfit, edm, errdef, nvpar, nparx, icstat);
  std::cout << "Best fit:" << std::endl;
  minuit->mnprin(4, lfit);
  delete fit;
  fit = nullptr;
  minuit = nullptr;

  // binary search through signal normalizations
  double r;
  double rcrit_low;
  double rcrit_high;
  float ns = 0;
  double converge = 0.5;
  double increment = 20.0;
  do {
    ns += increment;

    // conditional best fit
    Fit* fit = new Fit(signals, data);
    std::map<std::string, float> n = {{signal_name, ns}};
    std::map<std::string, bool> f = {{signal_name, true}};
    std::cout << "NS = " << ns << std::endl;
    minuit = (*fit)(e_range, r_range, &n, &f);
    double lfix;
    minuit->mnstat(lfix, edm, errdef, nvpar, nparx, icstat);
    if (VERBOSE) {
      std::cout << "Conditional best fit:" << std::endl;
      minuit->mnprin(4, lfix);
    }

    r = 2.0 * (lfix - lfit);  // r for this experiment
    if (VERBOSE) {
      std::cout << "r = " << r << std::endl;
    }

    // normalizations for fake data sets: sampled from ns and conditional best
    // fit backgrounds
    std::map<std::string, double> snormf;
    for (size_t i=0; i<signals.size(); i++) {
      double par, err;
      minuit->GetParameter(i, par, err);
      snormf[signals[i].name] = signals[i].rate * par;
    }
    snormf[signal_name] = ns;

    delete fit;
    fit = nullptr;
    minuit = nullptr;

    // generate fake experiments for this experiment to find C.L.% Rcrit
    std::vector<double> fake_rs;
    for (unsigned i=0; i<nmc; i++) {
      if (VERBOSE) {
        std::cout << "--------------------------------------------------------"
                  << std::endl;
      }
      std::cout << "Ns = " << ns << ", fake experiment " << i << std::endl;

      // fake data
      TNtuple* ntfake = gen.make_dataset(true, &snormf);
      Dataset fakedata("fakedata", ntfake);

      // best fit
      Fit* fit = new Fit(signals, fakedata);
      minuit = (*fit)(e_range, r_range);
      double lfit_fake;
      minuit->mnstat(lfit_fake, edm, errdef, nvpar, nparx, icstat);
      if (VERBOSE) {
        std::cout << "FD best fit:" << std::endl;
        minuit->mnprin(4, lfit_fake);
      }

      // conditional best fit
      n = {{signal_name, ns}};
      f = {{signal_name, true}};
      minuit = (*fit)(e_range, r_range, &n, &f);
      double lfix_fake;
      minuit->mnstat(lfix_fake, edm, errdef, nvpar, nparx, icstat);
      if (VERBOSE) {
        std::cout << "FD conditional best fit:" << std::endl;
        minuit->mnprin(4, lfix_fake);
      }

      delete fit;
      fit = nullptr;
      minuit = nullptr;

      double rfake = 2.0 * (lfix_fake - lfit_fake);
      std::cout << std::setprecision(9)
                << "lfit_fake = " << lfit_fake
                << ", lfix_fake = " << lfix_fake
                << ", rfake = " << rfake
                << ", ns = " << ns << std::endl;
      fake_rs.push_back(rfake);
      hr->Fill(ns, rfake);
    }

    // find R that C.L.% of Rs fall within
    std::sort(fake_rs.begin(), fake_rs.end());
    for (size_t k=0; k<fake_rs.size(); k++) {
      std::cout << fake_rs[k] << " ";
    }
    std::cout << std::endl;
    std::cout << std::floor(confidence * fake_rs.size()) << std::endl;
    rcrit_low = fake_rs[std::floor(1.0-confidence * fake_rs.size())];
    rcrit_high = fake_rs[std::floor(confidence * fake_rs.size())];

    std::cout << std::setprecision(9)
              << "lfit = " << lfit
              << ", lfix = " << lfix
              << ", r = " << r
              << ", rcrit low = " << rcrit_low
              << ", rcrit high = " << rcrit_high
              << ", ns = " << ns
              << std::endl;

    // experiment rejected, but increment is big: refine the binary search
    if (r > rcrit_high && increment > converge) {
      ns -= increment;
      increment /= 2;
    }

  } while (r < rcrit_high || increment > (2*converge));

  return ns;
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

  // find sensitivity limit
  TH2F hr("hr", "hr", 50, 0, 50, 10000, 150, 250);
  FakeDataGenerator gen(fc.signals, fc.e_range, fc.r_range);
  std::vector<double> rs;
  const unsigned nmc = fc.mc_trials;
  const unsigned nfake = fc.fake_experiments;
  for (unsigned i=0; i<nmc; i++) {
    float ul = get_ul(fc.signals, fc.signal_name, fc.confidence, gen, nfake,
                      &hr, fc.e_range, fc.r_range);
    rs.push_back(ul);
    std::cout << "Experiment " << i << ": UL N=" << ul << std::endl
              << "========================================================="
              << std::endl;
  }
  for (size_t i=0; i<rs.size(); i++) {
    std::cout << std::setprecision(9) << rs[i] << " " << std::endl;
  }
  std::cout << std::endl;

  // results
  double medrs = median(rs);
  std::cout << std::endl << "Median N_s: " << std::setprecision(8) << medrs
            << std::endl;

  TH1F* hsum = nullptr;
  TH1F* hsum_ns = nullptr;
  TH1F* hext = nullptr;
  TH1F* hcosmo = nullptr;

  // attempt to make a non-horrible plot
  gStyle->SetOptStat(0);
  std::vector<int> colors = { 1,  2,  3,  4,  6,  7,  8,  9, 11,
                              29, 5, 30, 34, 38, 40, 42, 45, 49 };
  std::map<std::string, std::vector<std::string> > chains = {
    {"External", {"av_bi214", "av_tl208", "water_bi214", "water_tl208", "pmt_bg"}},
    {"Cosmogenic", {"sb124", "sn126", "y88"}}
  };
  TCanvas c1;
  c1.SetLogy();
  TLegend stan(0.75, 0.55, 0.88, 0.88);
  stan.SetBorderSize(0);
  stan.SetFillColor(0);
  bool first = true;
  size_t icolor = 0;
  for (size_t i=0; i<fc.signals.size(); i++) {
    double p = (fc.signals[i].name == fc.signal_name ? medrs/fc.signals[i].rate : 1.0);
    fc.signals[i].histogram->Rebin(4);
    fc.signals[i].histogram->SetDirectory(0);
    fc.signals[i].histogram->GetXaxis()->SetRangeUser(2, 3.5);
    fc.signals[i].histogram->SetTitle("");
    fc.signals[i].histogram->SetXTitle("Energy (MeV)");
    char binsize[20];
    snprintf(binsize, 20, "%1.0f", fc.signals[i].histogram->GetXaxis()->GetBinWidth(1) * 1000);
    fc.signals[i].histogram->SetYTitle(("Counts per " + std::string(binsize) + " keV bin").c_str());
    fc.signals[i].histogram->Scale(fc.signals[i].rate * fc.live_time * p);
    fc.signals[i].histogram->SetAxisRange(1e-2, 1e2, "Y");

    std::vector<std::string> ve = chains["External"];
    std::vector<std::string> vc = chains["Cosmogenic"];
    if (std::find(ve.begin(), ve.end(), fc.signals[i].name) != ve.end()) {
      if (hext == nullptr) {
        hext = (TH1F*) fc.signals[i].histogram->Clone("hext");
        hext->SetLineColor(colors[icolor]);
        hext->SetLineWidth(2);
        stan.AddEntry(hext, "External");
        icolor++;
      }
      hext->Add(fc.signals[i].histogram);
    }
    else if (std::find(vc.begin(), vc.end(), fc.signals[i].name) != vc.end()) {
      if (hcosmo == nullptr) {
        hcosmo = (TH1F*) fc.signals[i].histogram->Clone("hcosmo");
        hcosmo->SetLineColor(colors[icolor]);
        hcosmo->SetLineWidth(2);
        stan.AddEntry(hcosmo, "Cosmogenic");
        icolor++;
      }
      hcosmo->Add(fc.signals[i].histogram);
    }
    else {
      if (hsum == nullptr) {
        hsum = (TH1F*) fc.signals[i].histogram->Clone("hsum");
        hsum->SetLineColor(1);
        hsum->SetLineWidth(2);
        stan.AddEntry(hsum, "Sum");
      }
      if (hsum_ns == nullptr) {
        hsum_ns = (TH1F*) fc.signals[i].histogram->Clone("hsum_ns");
        hsum_ns->SetLineColor(1);
        hsum_ns->SetLineStyle(2);
        hsum_ns->SetLineWidth(2);
        stan.AddEntry(hsum_ns, "Sum, 0 signal");
      }

      fc.signals[i].histogram->SetLineColor(colors[icolor]);
      fc.signals[i].histogram->SetLineWidth(2);
      fc.signals[i].histogram->Draw(first ? "": "same");
      stan.AddEntry(fc.signals[i].histogram, fc.signals[i].title.c_str());

      icolor++;
      first = false;
    }

    hsum->Add(fc.signals[i].histogram);
    if (fc.signals[i].name != fc.signal_name) {
      hsum_ns->Add(fc.signals[i].histogram);
    }
  }

  hext->Draw("same");
  hcosmo->Draw("same");
  hsum_ns->Draw("same");
  hsum->Draw("same");
  
  stan.Draw();
  c1.Update();
  c1.SaveAs("fit_spectrum_a.C");
  c1.SaveAs("fit_spectrum_a.pdf");

  TCanvas c2;
  c2.cd();
  hr.Draw("col z");
  c2.SaveAs("hr.C");
  c2.SaveAs("hr.pdf");

  return 0;
}

