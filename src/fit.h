#ifndef __FIT_H__
#define __FIT_H__

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <TChain.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>
#include "ll.h"

class TMinuit;


// a concatenation of ntuples
class Dataset {
  public:
    // construct using ntuple files on disk
    Dataset(std::string _name, std::vector<std::string> filenames) : name(_name) {
      TChain* _events = new TChain("tev");
      this->events = dynamic_cast<TTree*>(_events);

      for (size_t i=0; i<filenames.size(); i++) {
        std::cout << "Adding: " << filenames[i] << std::endl;
        _events->Add(filenames[i].c_str());
      }

      this->events->SetBranchStatus("*", 0);
      this->events->SetBranchStatus("e", 1);
      this->events->SetBranchStatus("r", 1);
      this->nevents = events->GetEntries();
      std::cout << "Dataset " << this->name << ": " << this->nevents << " events" << std::endl;
    }

    // construct using a ttree in memory
    Dataset(std::string _name, TTree* _events) : name(_name), events(_events) {
      this->events->SetBranchStatus("*", 0);
      this->events->SetBranchStatus("e", 1);
      this->events->SetBranchStatus("r", 1);
      this->nevents = events->GetEntries();
      std::cout << "Dataset " << this->name << ": " << this->nevents << " events" << std::endl;
    }

    virtual ~Dataset() {
      delete this->events;
    }

    std::string name;
    size_t nevents;
    TTree* events;
};


// a structure for passing information to minuit
struct Experiment {
  size_t nevents;  // total event count in dataset (rows in lut)
  const std::vector<Signal>* signals;
  GPULL* ll;
};


// user interface to the MINUIT ML fit
class Fit {
  public:
    Fit(const std::vector<Signal>& signals, const Dataset& data);

    virtual ~Fit() {
      delete this->ll;
      delete this->lut;
      delete this->minuit;
    }

    // run the fit, optionally specifying some different starting normalizations
    // and fixed/unfixed status
    TMinuit* operator()(std::map<std::string, float>* _norms=nullptr, std::map<std::string, bool>* _fix=nullptr);

  protected:
    // construct a P(x) lookup table
    float* build_lut(const std::vector<Signal>& signals, const Dataset& data) const;

    // minuit fit function
    static void nll(int& ndim, double* gout, double& result, double* par, int flags);

    // global container for passing state to the minuit fit function
    static Experiment experiment;

  private:
    float* lut;  // P(x) lookup table
    GPULL* ll;  // GPU LL interface
    TMinuit* minuit;
    std::vector<Signal> signals;
};

#endif  // __FIT_H__

