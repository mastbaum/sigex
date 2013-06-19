#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <assert.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TMinuit.h>

#include "util.h"
#include "fit.h"

//#define DEBUG  // enable debugging printout

#ifndef VERBOSE
#define VERBOSE false
#endif

std::vector<Signal> Fit::signals;
Range<float> Fit::e_range;
Range<float> Fit::r_range;
float* Fit::norms;
TH1* Fit::data;


Fit::Fit(const std::vector<Signal>& signals, const Dataset& data) {
  Fit::signals = signals;

  // bin the data
  Fit::data = (TH1*) signals[0].histogram->Clone("data");
  Fit::data->Reset();
  if (Fit::data->IsA() == TH2F::Class()) {
    data.events->Draw("r:e>>data");
  }
  else {
    data.events->Draw("e>>data");
  }
  std::cout << "data: " << Fit::data->Integral() << std::endl;

  // set up minuit instance
  minuit = new TMinuit(signals.size());
  int eflag;
  double p1 = -1;
  minuit->mnexcm("SET PRINTOUT", &p1, 1, eflag);
  p1 = 0.5;
  minuit->mnexcm("SET ERR", &p1, 1, eflag);
  p1 = 2;
  minuit->mnexcm("SET STR", &p1, 1, eflag);
  
  minuit->SetFCN(Fit::nll);

  for (size_t i=0; i<signals.size(); i++) {
    this->minuit->mnparm(i, signals[i].name.c_str(), 1.0, 1e-3, 0, 0, eflag);
  }
}


void Fit::nll(int& ndim, double* gout, double& result, double* par, int flags) {
  size_t n = Fit::signals.size();
  result = 0;

  // make summed fit histogram
  TH1* hfit = (TH1*) signals[0].histogram->Clone("hfit");
  hfit->Reset();
  for (size_t i=0; i<n; i++) {
    TH1* h = Fit::signals.at(i).histogram;
    h->Scale(1.0/h->Integral());
    hfit->Add(h, Fit::norms[i] * par[i]);
  }

  // loop over bins
  if (hfit->IsA() == TH2F::Class()) {
    TH2F* h2 = dynamic_cast<TH2F*>(hfit);
    for (int i=1; i<h2->GetNbinsX(); i++) {
      if (i < h2->GetXaxis()->FindBin(Fit::r_range.min) || i > h2->GetXaxis()->FindBin(Fit::r_range.max)) {
        continue;
      }
      for (int j=1; j<h2->GetNbinsY(); j++) {
        if (j < h2->GetYaxis()->FindBin(Fit::e_range.min) || j > h2->GetYaxis()->FindBin(Fit::e_range.max)) {
          continue;
        }
        double nexp = h2->GetBinContent(i, j);
        double nobs = dynamic_cast<TH2F*>(Fit::data)->GetBinContent(i, j);
        result += (nexp - nobs * TMath::Log(TMath::Max(1e-12, nexp)));
        //std::cout << "- " << nexp << " " << nobs << " " << TMath::Log(nexp) << " // " << result << std::endl;
      }
    }
  }
  else {
    for (int i=1; i<hfit->GetNbinsX(); i++) {
      if (i < hfit->FindBin(Fit::e_range.min) || i > hfit->FindBin(Fit::e_range.max)) {
        continue;
      }
      double nexp = hfit->GetBinContent(i);
      double nobs = Fit::data->GetBinContent(i);
      result += (nexp - nobs * TMath::Log(TMath::Max(1e-12, nexp)));
      //std::cout << "- " << nexp << " " << nobs << " " << TMath::Log(nexp) << " // " << result << std::endl;
    }
  }

  // constraints
  for (size_t i=0; i<n; i++) {
    if (Fit::signals.at(i).constraint > 0) {
      result += 0.5 * TMath::Power((par[i]-1.0)
                / Fit::signals.at(i).constraint, 2);
    }
  }
  delete hfit;

#ifdef DEBUG
  // print parameters at each iteration
  std::cout << "+ ";
  for (size_t i=0; i<n; i++) {
    std::cout << par[i] * Fit::signals.at(i).rate << " (" << par[i] << ") \t";
  }
  std::cout << result << std::endl;
#endif
}


TMinuit* Fit::operator()(Range<float> _e_range, Range<float> _r_range,
                         float live_time,
                         std::map<std::string, float>* _norms,
                         std::map<std::string, bool>* _fix,
                         bool run_minos) {
  Fit::e_range = _e_range;
  Fit::r_range = _r_range;

  // build a complete list of normalizations
  float* norms = new float[this->signals.size()];
  for (size_t i=0; i<this->signals.size(); i++) {
    norms[i] = live_time * this->signals[i].rate;
    if (_norms != nullptr) {
      auto new_norm = _norms->find(this->signals[i].name);
      if (new_norm != _norms->end()) {
        norms[i] = live_time * new_norm->second;
        this->signals[i].rate = live_time * new_norm->second;
      }
    }
  }
  delete Fit::norms;
  Fit::norms = norms;

  // run fit
  if (VERBOSE) {
    std::cout << "Starting fit..." << std::endl;
  }
  TStopwatch timer;
  timer.Start();

  // set parameters fixed
  int eflag;
  for (size_t i=0; i<this->signals.size(); i++) {
    if (signals[i].fixed) {
      this->minuit->FixParameter(i);
    }
    else {
      this->minuit->Release(i);
    }
    if (_fix != nullptr) {
      auto new_fix = _fix->find(this->signals[i].name);
      if (new_fix != _fix->end()) {
        if (new_fix->second) {
          double arglist[2];
          int ierr = 0;
          arglist[0] = i + 1;
          arglist[1] = 1.0;
          this->minuit->mnexcm("SET PAR", arglist, 2, ierr);
          this->minuit->FixParameter(i);
        }
        else {
          this->minuit->Release(i);
        }
      }
    }
  }

  minuit->mnexcm("SIMPLEX 4000", 0, 0, eflag);
  minuit->mnexcm("MIGRAD 4000", 0, 0, eflag);
  if (run_minos) {
    minuit->mnexcm("MINOS", 0, 0, eflag);
  }

  if (VERBOSE) {
    std::cout << "Fit completed in " << timer.RealTime() << " seconds"
              << std::endl;

    double sum = 0;
    for (size_t i=0; i<this->signals.size(); i++) {
      double p, err;
      minuit->GetParameter(i, p, err);
      std::cout << signals[i].name << ": " << p << " * " << norms[i] << " = " << p * norms[i] << std::endl;
      sum += p * norms[i];
    }
    std::cout << "total fit events: " << sum << std::endl;
  }

  return minuit;
}

