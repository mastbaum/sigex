#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <assert.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TMinuit.h>
#include "util.h"
#include "fit.h"

// create fake datasets by sampling histograms
class FakeDataGenerator {
  public:
    FakeDataGenerator(std::vector<Signal> signals, Range<float> _e_range, Range<float> _r_range)
      : e_range(_e_range), r_range(_r_range) {
      for (auto it=signals.cbegin(); it!=signals.cend(); it++) {
        this->pdfs[it->name] = it->histogram;
        this->default_norms[it->name] = it->rate;
      }
    }
    virtual ~FakeDataGenerator() {}

    // create a dataset: an ntuple with fields "r:e"
    // if histograms are TH2, both are filled; if TH1, r is set to 0.
    TNtuple* make_dataset(bool poisson=true, std::map<std::string, double>* _norms=nullptr);

  protected:
    std::map<std::string, TH1*> pdfs;
    std::map<std::string, double> default_norms;
    Range<float> e_range;
    Range<float> r_range;
};


TNtuple* FakeDataGenerator::make_dataset(bool poisson, std::map<std::string, double>* _norms) {
  std::cout << "FakeDataGenerator::make_dataset: Generating dataset..." << std::endl;
  TNtuple* nt = new TNtuple("tev", "events", "r:e");
  double r = 0;
  double e = 0;
  for (auto it=this->pdfs.begin(); it!=this->pdfs.end(); it++) {
    int nexpected = ((_norms && _norms->find(it->first) != _norms->end()) ? (*_norms)[it->first] : this->default_norms[it->first]);

    double integral;
    if (it->second->IsA() == TH2F::Class()) {
      TH2F* ht = dynamic_cast<TH2F*>(it->second);
      int x1 = it->second->GetXaxis()->FindBin(r_range.min);
      int x2 = it->second->GetXaxis()->FindBin(r_range.max);
      int y1 = it->second->GetYaxis()->FindBin(e_range.min);
      int y2 = it->second->GetYaxis()->FindBin(e_range.max);
      integral = ht->Integral(x1, x2, y1, y2);
      if (integral == 0) {
        continue;
      }
      int nobserved = nexpected;
      nexpected /= integral;
      if (poisson && nexpected > 0) {
        nobserved = gRandom->Poisson(nexpected);
      }
      for (int i=0; i<nobserved; i++) {
        do {
          ht->GetRandom2(r, e);
        } while(e > e_range.max || e < e_range.min || r > r_range.max || r < r_range.min);
        nt->Fill(r, e);
      }
      std::cout << "FakeDataGenerator::make_dataset: " << it->first << ": " << nobserved << " events" << std::endl;
    }
    else if (it->second->IsA() == TH1D::Class()) {
      TH1D* ht = dynamic_cast<TH1D*>(it->second);
      int x1 = it->second->GetXaxis()->FindBin(e_range.min);
      int x2 = it->second->GetXaxis()->FindBin(e_range.max);
      integral = ht->Integral(x1, x2);
      if (integral == 0) {
        continue;
      }
      int nobserved = nexpected;
      nexpected /= integral;
      if (poisson && nexpected > 0) {
        nobserved = gRandom->Poisson(nexpected);
      }
      for (int i=0; i<nobserved; i++) {
        do {
          e = ht->GetRandom();
        } while(e > e_range.max || e < e_range.min);
        nt->Fill(r, e);
      }
      std::cout << "FakeDataGenerator::make_dataset: " << it->first << ": " << nobserved << " events" << std::endl;
    }
    else {
      std::cerr << "FakeDataGenerator::make_dataset: Unknown histogram class " << it->second->ClassName() << std::endl;
      assert(false);
    }
  }

  return nt;
}


// run an experiment, get an upper limit
float get_ul(std::vector<Signal> _signals, std::string signal_name, float confidence, FakeDataGenerator& gen, unsigned nmc) {
  std::vector<Signal> signals = _signals;

  // make experiment data
  std::map<std::string, double> snorm = {{signal_name, 0}};
  TNtuple* nt = gen.make_dataset(true, &snorm);
  Dataset data("data", nt);

  //TH1F* hdata = new TH1F("hdata", "hdata", 300, 2, 5);
  //nt->Draw("e>>hdata");

  // best fit
  Fit* fit = new Fit(signals, data);
  TMinuit* minuit = (*fit)();
  double lfit, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat(lfit, edm, errdef, nvpar, nparx, icstat);
  std::cout << "Best fit:" << std::endl;
  minuit->mnprin(4, lfit);
  //TH1F* hbf = histogram(signals, minuit, "hbf", 300, 2, 5);
  delete fit;
  fit = nullptr;
  minuit = nullptr;

  // binary search through signal normalizations
  double r;
  double rcrit_low;
  double rcrit_high;
  float ns = 0;
  double converge = 2; //0.5;
  double increment = 20.0;
  do {
    ns += increment;

    // conditional best fit
    Fit* fit = new Fit(signals, data);
    std::map<std::string, float> n = {{signal_name, ns}};
    std::map<std::string, bool> f = {{signal_name, true}};
    std::cout << "NS = " << ns << std::endl;
    minuit = (*fit)(&n, &f);
    double lfix;
    minuit->mnstat(lfix, edm, errdef, nvpar, nparx, icstat);
    std::cout << "Conditional best fit:" << std::endl;
    minuit->mnprin(4, lfix);

    r = 2.0 * (lfix - lfit);  // r for this experiment
    std::cout << "r = " << r << std::endl;

    // normalizations for fake data sets: sampled from ns and conditional best fit backgrounds
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
      std::cout << "---------------------------------------------------------" << std::endl;
      std::cout << "Ns = " << ns << ", fake experiment " << i << std::endl;

      // fake data
      TNtuple* ntfake = gen.make_dataset(true, &snormf);
      Dataset fakedata("fakedata", ntfake);

      // best fit
      Fit* fit = new Fit(signals, fakedata);
      minuit = (*fit)();
      double lfit_fake;
      minuit->mnstat(lfit_fake, edm, errdef, nvpar, nparx, icstat);
      std::cout << "FD best fit:" << std::endl;
      minuit->mnprin(4, lfit_fake);

      // conditional best fit
      n = {{signal_name, ns}};
      f = {{signal_name, true}};
      minuit = (*fit)(&n, &f);
      double lfix_fake;
      minuit->mnstat(lfix_fake, edm, errdef, nvpar, nparx, icstat);
      std::cout << "FD conditional best fit:" << std::endl;
      minuit->mnprin(4, lfix_fake);

      delete fit;
      fit = nullptr;
      minuit = nullptr;

      double rfake = 2.0 * (lfix_fake - lfit_fake);
      std::cout << std::setprecision(9) << "lfit_fake = " << lfit_fake << ", lfix_fake = " << lfix_fake << ", rfake = " << rfake << ", ns = " << ns << std::endl;
      fake_rs.push_back(rfake);
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

    std::cout << std::setprecision(9) << "lfit = " << lfit << ", lfix = " << lfix << ", r = " << r << ", rcrit low = " << rcrit_low << ", rcrit high = " << rcrit_high << ", ns = " << ns << std::endl;

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
  FakeDataGenerator gen(fc.signals, fc.e_range, fc.r_range);
  std::vector<double> rs;
  const unsigned nmc = fc.mc_trials;
  for (unsigned i=0; i<nmc; i++) {
    float ul = get_ul(fc.signals, fc.signal_name, fc.confidence, gen, nmc);
    rs.push_back(ul);
    std::cout << "Experiment " << i << ": UL N=" << ul << std::endl;
    std::cout << "=========================================================" << std::endl;
  }
  for (size_t i=0; i<rs.size(); i++) {
    std::cout << std::setprecision(9) << rs[i] << " " << std::endl;
  }
  std::cout << std::endl;

  // results
  double medrs = median(rs);
  std::cout << std::endl << "Median N_s: " << std::setprecision(8) << medrs << std::endl;

  TH1F* hsum = NULL;
  gStyle->SetOptStat(0);
  std::vector<int> colors = {1, 2, 3, 4, 6, 7, 8, 9, 11, 29, 5, 30, 34, 38, 40, 42, 45, 49};
  TCanvas c1;
  c1.SetLogy();
  TLegend stan(0.75, 0.55, 0.88, 0.88);
  stan.SetBorderSize(0);
  stan.SetFillColor(0);
  bool first = true;
  for (size_t i=0; i<fc.signals.size(); i++) {
    double p = (fc.signals[i].name == fc.signal_name ? medrs/fc.signals[i].rate : 1.0);
    fc.signals[i].histogram->SetLineColor(colors[i]);
    fc.signals[i].histogram->SetLineWidth(2);
    fc.signals[i].histogram->Rebin(4);
    fc.signals[i].histogram->SetDirectory(0);
    fc.signals[i].histogram->GetYaxis()->SetRangeUser(1e-3, 1e2);
    fc.signals[i].histogram->GetXaxis()->SetRangeUser(2, 3.5);
    fc.signals[i].histogram->SetTitle("");
    fc.signals[i].histogram->SetXTitle("Energy (MeV)");
    char binsize[20];
    snprintf(binsize, 20, "%1.0f", fc.signals[i].histogram->GetXaxis()->GetBinWidth(1) * 1000);
    fc.signals[i].histogram->SetYTitle(("Counts per " + std::string(binsize) + " keV bin").c_str());
    fc.signals[i].histogram->Scale(fc.signals[i].rate * fc.live_time * p);
    fc.signals[i].histogram->Draw(first ? "" : "same");
    stan.AddEntry(fc.signals[i].histogram, fc.signals[i].title.c_str());
    std::cout << "Adding " << fc.signals[i].name << ", " << fc.signals[i].histogram->GetNbinsX() << " bins)" << std::endl;
    if (first) {
      hsum = (TH1F*)fc.signals[i].histogram->Clone("hsum");
      hsum->Reset();
      hsum->SetLineWidth(2);
      hsum->SetMarkerStyle(20);
      std::cout << "hsum has " << hsum->GetNbinsX() << " bins" << std::endl;
    }
    hsum->Add(fc.signals[i].histogram);
    first = false;
  }
  hsum->Draw("same");
  
  stan.AddEntry(hsum, "MC");
  stan.Draw();
  c1.Update();
  c1.SaveAs("fit_spectrum.C");
  c1.SaveAs("fit_spectrum.pdf");

  return 0;
}

