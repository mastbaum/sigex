#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <assert.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TMinuit.h>
#include "util.h"
#include "ll.h"
#include "fit.h"


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

    TNtuple* make_dataset(bool poisson=true, std::map<std::string, double>* _norms=nullptr) {
      std::cout << "FakeDataGenerator::make_dataset: Generating dataset..." << std::endl;
      TNtuple* nt = new TNtuple("tev", "events", "r:e");
      double r = 0;
      double e = 0;
      for (auto it=this->pdfs.begin(); it!=this->pdfs.end(); it++) {
        int nexpected = (_norms && _norms->find(it->first) != _norms->end() ? _norms->at(it->first) : this->default_norms[it->first]);
        int nobserved = nexpected;
        if (poisson && nexpected > 0) {
          nobserved = gRandom->Poisson(nexpected);
        }
        std::cout << "FakeDataGenerator::make_dataset: " << it->first << ": " << nobserved << " events" << std::endl;

        int x1 = it->second->GetXaxis()->FindBin(r_range.min);
        int x2 = it->second->GetXaxis()->FindBin(r_range.max);
        int y1 = it->second->GetYaxis()->FindBin(e_range.min);
        int y2 = it->second->GetYaxis()->FindBin(e_range.max);

        double integral;
        if (it->second->IsA() == TH2F::Class()) {
          TH2F* ht = dynamic_cast<TH2F*>(it->second);
          integral = ht->Integral(x1, x2, y1, y2);
          if (integral == 0) {
            continue;
          }
          for (int i=0; i<nobserved; i++) {
            do {
              ht->GetRandom2(r, e);
            } while(e > e_range.max || e < e_range.min || r > r_range.max || r < r_range.min);
            nt->Fill(r, e);
          }
        }
        else if (it->second->IsA() == TH1D::Class()) {
          TH1D* ht = dynamic_cast<TH1D*>(it->second);
          integral = ht->Integral(y1, y2);
          if (integral == 0) {
            continue;
          }
          for (int i=0; i<nobserved; i++) {
            do {
              e = ht->GetRandom();
            } while(e > e_range.max || e < e_range.min);
            nt->Fill(r, e);
          }
        }
        else {
          std::cerr << "FakeDataGenerator::make_dataset: Unknown histogram class " << it->second->ClassName() << std::endl;
          assert(false);
        }
      }
      return nt;
    }

  protected:
    std::map<std::string, TH1*> pdfs;
    std::map<std::string, double> default_norms;
    Range<float> e_range;
    Range<float> r_range;
};


// global for passing configuration to minuit
Experiment experiment;


float* build_lut(const std::vector<Signal>& signals, const Dataset& data) {
  int nevents = data.nevents;
  float* lut = new float[signals.size() * nevents];

  std::vector<float> minima;
  for (auto it=signals.begin(); it!=signals.end(); ++it) {
    minima.push_back((*it).histogram->GetMinimum(0) * 0.001);
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
        std::cerr << "build_lut: Unknown histogram class " << signals[j].histogram->ClassName() << std::endl;
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


void nll_gpu(int& ndim, double* gout, double& result, double* par, int flags) {
  size_t n = experiment.signals->size();
  result = 0;

  // sum(N) + constraints
  for (size_t i=0; i<n; i++) {
    result += par[i] * experiment.signals->at(i).rate;
    if (experiment.signals->at(i).constraint > 0) {
      result += 0.5 * TMath::Power((par[i]-1.0)/experiment.signals->at(i).constraint, 2);
    }
  }

  // sum(log(sum(N_j * P_j(x_i))))
  result -= (*experiment.ll)(par);

#ifdef DEBUG
  // print parameters at each iteration
  std::cout << "+ ";
  for (size_t i=0; i<n; i++) {
    std::cout << par[i] * experiment.signals->at(i).norm << " (" << par[i] << ") " << "\t";
  }
  std::cout << result << std::endl;
#endif
}


// if _lut is null, compute lut and set _lut to point to it. if not null, use it as the lut.
// for reuse in subsequent fits. be sure to delete it.
TMinuit* run_fit(std::vector<Signal>& signals, const Dataset& data, float*& _lut) {
  float* lut;
  if (_lut == nullptr) {
    std::cout << "Building LUT..." << std::endl;
    lut = build_lut(signals, data);
    _lut = lut;
  }
  else {
    lut = _lut;
  }

  float* norms = new float[signals.size()];
  for (size_t i=0; i<signals.size(); i++) {
    norms[i] = signals[i].rate;
    //std::cout << signals[i].name << ": " << norms[i] << std::endl;
  }

  GPULL* ll = new GPULL(norms, lut, signals.size(), data.nevents);
  std::cout << "nevents = " << data.nevents << std::endl;

  // set globals for minuit
  experiment.ll = ll;
  experiment.nevents = data.nevents;
  experiment.signals = &signals;

  // run fit
  std::cout << "Starting fit..." << std::endl;
  TStopwatch timer;
  timer.Start();

  TMinuit* minuit = new TMinuit(signals.size());
  int eflag;
  double p1 = -1;
  minuit->mnexcm("SET PRINTOUT", &p1, 1, eflag);
  p1 = 0.5;
  minuit->mnexcm("SET ERR", &p1, 1, eflag);
  p1 = 2;
  minuit->mnexcm("SET STR", &p1, 1, eflag);
  
  minuit->SetFCN(nll_gpu);

  for (size_t i=0; i<signals.size(); i++) {
    minuit->mnparm(i, signals[i].name.c_str(), 1.0, 0.01, 0, 100, eflag);
    if (signals[i].fixed) {
      minuit->FixParameter(i);
    }
  }

  minuit->mnexcm("SIMPLEX", 0, 0, eflag);
  minuit->mnexcm("MIGRAD", 0, 0, eflag);
  minuit->mnexcm("MINOS", 0, 0, eflag);

  delete ll;
  delete[] norms;

  std::cout << "Fit completed in " << timer.RealTime() << " seconds" << std::endl;

  return minuit;
}


// run an experiment, get a 90% upper limit
double get_ul(std::vector<Signal> _signals, std::string signal_name, float confidence, FakeDataGenerator& gen, unsigned nmc) {
  std::vector<Signal> signals = _signals;

  // make experiment data
  std::map<std::string, double> snorm = {{signal_name, 0}};
  TNtuple* nt = gen.make_dataset(true, &snorm);
  Dataset data("data", nt);

  // best fit
  float* lut = nullptr;
  TMinuit* minuit = run_fit(signals, data, lut);
  double lfit, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat(lfit, edm, errdef, nvpar, nparx, icstat);
  std::cout << "Best fit:" << std::endl;
  minuit->mnprin(4, lfit);

  /*
  double signal_bestfit = 0;
  double err = 0;
  for (size_t i=0; i<signals.size(); i++) {
    if (signals[i].name == signal_name) {
      minuit->GetParameter(i, signal_bestfit, err);
    }
  }
  */

  delete minuit;
  minuit = nullptr;

  delete lut;
  lut = nullptr;

  // binary search through signal normalizations
  double r;
  double rcrit;
  double ns = 0;
  double converge = 0.5;
  double increment = 100.0;
  do {
    ns += increment;

    // set and fix signal normalization
    for (size_t i=0; i<signals.size(); i++) {
      if (signals[i].name == signal_name) {
        signals[i].fixed = true;
        signals[i].rate = ns;
      }
    }

    // conditional best fit
    minuit = run_fit(signals, data, lut);
    double lfix;
    minuit->mnstat(lfix, edm, errdef, nvpar, nparx, icstat);
    std::cout << "Conditional best fit:" << std::endl;
    minuit->mnprin(4, lfit);

    r = 2.0 * (lfix - lfit);  // r for this experiment
    std::cout << "r = " << r << std::endl;

    // normalizations for fake data sets, sampled from ns and conditional best fit backgrounds
    std::map<std::string, double> snormf;
    for (size_t i=0; i<signals.size(); i++) {
      double par, err;
      minuit->GetParameter(i, par, err);
      snormf[signals[i].name] = signals[i].rate * par;
    }
    snormf[signal_name] = ns;

    delete minuit;
    minuit = nullptr;

    delete lut;
    lut = nullptr;

    // generate fake experiments for this experiment to find 90% Rcrit
    std::vector<double> fake_rs;
    for (unsigned i=0; i<nmc; i++) {
      std::cout << "Ns = " << ns << ", fake experiment " << i << std::endl;

      // fake data
      TNtuple* ntfake = gen.make_dataset(true, &snormf);
      Dataset fakedata("fakedata", ntfake);

      // best fit
      for (size_t i=0; i<signals.size(); i++) {
        if (signals[i].name == signal_name) {
          signals[i].fixed = false;
        }
      }

      minuit = run_fit(signals, fakedata, lut);
      double lfit_fake;
      minuit->mnstat(lfit_fake, edm, errdef, nvpar, nparx, icstat);
      std::cout << "FD best fit:" << std::endl;
      minuit->mnprin(4, lfit_fake);

      delete minuit;
      minuit = nullptr;

      // conditional best fit
      for (size_t i=0; i<signals.size(); i++) {
        if (signals[i].name == signal_name) {
          signals[i].fixed = true;
        }
      }

      minuit = run_fit(signals, fakedata, lut);
      double lfix_fake;
      minuit->mnstat(lfix_fake, edm, errdef, nvpar, nparx, icstat);
      std::cout << "FD conditional best fit:" << std::endl;
      minuit->mnprin(4, lfix_fake);

      delete minuit;
      minuit = nullptr;

      delete lut;
      lut = nullptr;

      double rfake = 2.0 * (lfix_fake - lfit_fake);
      std::cout << std::setprecision(9) << "lfit_fake = " << lfit_fake << ", lfix_fake = " << lfix_fake << ", rfake = " << rfake << ", ns = " << ns << std::endl;
      if (rfake <= 0) {
        continue;
      }
      fake_rs.push_back(rfake);
    }

    // find R that C.L.% of Rs fall below
    std::sort(fake_rs.begin(), fake_rs.end());
    for (size_t k=0; k<fake_rs.size(); k++) {
      std::cout << fake_rs[k] << " ";
    }
    std::cout << std::endl;
    std::cout << std::floor(confidence * fake_rs.size()) << std::endl;
    rcrit = fake_rs[std::floor(confidence * fake_rs.size())];

    std::cout << std::setprecision(9) << "lfit = " << lfit << ", lfix = " << lfix << ", r = " << r << ", rcrit = " << rcrit << ", ns = " << ns << std::endl;

    // experiment rejected, but increment is big: refine the binary search
    if (r > rcrit && increment > converge) {
      ns -= increment;
      increment /= 2;
    }

  } while (r < rcrit || increment > (2*converge));

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

  // sensitivity
  FakeDataGenerator gen(fc.signals, fc.e_range, fc.r_range);
  std::vector<double> rs;
  const unsigned nmc = fc.mc_trials;
  for (unsigned i=0; i<nmc; i++) {
    double ul = get_ul(fc.signals, fc.signal_name, fc.confidence, gen, nmc);
    std::cout << "Experiment " << i << ": UL N=" << ul << std::endl;
    rs.push_back(ul);
  }
  for (size_t i=0; i<rs.size(); i++) {
    std::cout << std::setprecision(9) << rs[i] << " " << std::endl;
  }
  std::cout << std::endl;

  double medrs = median(rs);
  std::cout << "Median N_s: " << std::setprecision(8) << medrs << std::endl;

  return 0;
}

