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
#include "ll.h"
#include "fit.h"

#define __GPU__  // use the GPU to accelerate NLL calculation
//#define DEBUG  // enable debugging printout

#ifndef VERBOSE
#define VERBOSE false
#endif

std::vector<Signal> Fit::signals;
GPULL* Fit::ll;
size_t Fit::nevents;


Fit::Fit(const std::vector<Signal>& signals, const Dataset& data) {
  Fit::signals = signals;
  Fit::nevents = data.nevents;

  this->lut = this->build_lut(Fit::signals, data);
  Fit::ll = new GPULL(this->lut, Fit::signals.size(), Fit::nevents);

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
    this->minuit->mnparm(i, signals[i].name.c_str(), 1.0, 1e-3, -50, 50, eflag);
  }
}


float* Fit::build_lut(const std::vector<Signal>& signals, const Dataset& data) {
  int nevents = data.nevents;
  float* lut = new float[signals.size() * nevents];

  std::vector<float> minima;
  for (auto it=signals.begin(); it!=signals.end(); ++it) {
    minima.push_back((*it).histogram->GetMinimum(0) * 0.0001);
  }

  float e;
  float r;
  data.events->SetBranchAddress("e", &e);
  data.events->SetBranchAddress("r", &r);

  for (int i=0; i<nevents; i++) {
    data.events->GetEntry(i);
    for (size_t j=0; j<signals.size(); j++) {
      float v = 0;
      if (signals[j].histogram->IsA() == TH2F::Class()) {
        v = dynamic_cast<TH2F*>(signals[j].histogram)->Interpolate(r, e);
      }
      else if (signals[j].histogram->IsA() == TH1D::Class()) {
        v = dynamic_cast<TH1D*>(signals[j].histogram)->Interpolate(e);
      }
      else {
        std::cerr << "build_lut: Unknown histogram class "
                  << signals[j].histogram->ClassName() << std::endl;
        assert(false);
      }
      
      if (v <= 0) {
        v = minima[j];
      }
      lut[i*signals.size() + j] = v;
    }
  }

  return lut;
}


void Fit::nll(int& ndim, double* gout, double& result, double* par,
              int flags) {
  size_t n = Fit::signals.size();
  result = 0;

  // sum(N) + constraints
  for (size_t i=0; i<n; i++) {
    result += par[i] * Fit::signals.at(i).rate;
    if (Fit::signals.at(i).constraint > 0) {
      result += 0.5 * TMath::Power((par[i]-1.0)
                / Fit::signals.at(i).constraint, 2);
    }
  }

  // fixme just delete this
  // -log(N!)
  int nevents = Fit::nevents;
  if (nevents < 200) {
    result -= TMath::Log(TMath::Factorial(nevents));
  }
  else {
    result -= TMath::Log(nevents * TMath::Log(nevents) - nevents);  // stirling
  }

  // sum(log(sum(N_j * P_j(x_i))))
#ifdef __GPU__
  result -= (*Fit::ll)(par);
#else
  std::cerr << "CPU NLL not yet implemented" << std::endl;
  assert(0);
#endif

#ifdef DEBUG
  // print parameters at each iteration
  std::cout << "+ ";
  for (size_t i=0; i<n; i++) {
    std::cout << par[i] * Fit::signals->at(i).norm << " (" << par[i] << ") \t";
  }
  std::cout << result << std::endl;
#endif
}


TMinuit* Fit::operator()(Range<float> e_range,
                         std::map<std::string, float>* _norms,
                         std::map<std::string, bool>* _fix,
                         bool run_minos) {
  // build a complete list of normalizations
  float* norms = new float[this->signals.size()];
  for (size_t i=0; i<this->signals.size(); i++) {
    int x1 = this->signals[i].histogram->FindBin(e_range.min);
    int x2 = this->signals[i].histogram->FindBin(e_range.max);
    double ii = this->signals[i].histogram->Integral(x1, x2);
    norms[i] = this->signals[i].rate / (ii / this->signals[i].histogram->Integral());
    if (_norms != nullptr) {
      auto new_norm = _norms->find(this->signals[i].name);
      if (new_norm != _norms->end()) {
        norms[i] = new_norm->second;
        this->signals[i].rate = new_norm->second;
      }
    }
  }

  // copy new normalizations to gpu
  this->ll->set_norms(norms, this->signals.size());

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

  // run fit
  if (VERBOSE) {
    std::cout << "Starting fit..." << std::endl;
  }
  TStopwatch timer;
  timer.Start();

  minuit->mnexcm("SIMPLEX 4000", 0, 0, eflag);
  minuit->mnexcm("MINIMIZE 4000", 0, 0, eflag);
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
      std::cout << signals[i].name << ": " << p * norms[i] << std::endl;
      sum += p * norms[i];
    }
    std::cout << "total fit events: " << sum << std::endl;
  }

  delete[] norms;
  norms = nullptr;

  return minuit;
}

