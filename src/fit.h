#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <TChain.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>

class TMinuit;
class GPULL;

// a concatenation of ntuples
class Dataset {
  public:
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


// construct an events * pdfs likelihood lookup table
// this is stored as a contiguous 1d array indexed by "iev*nsignals+isignal"
float* build_lut(const std::vector<Signal>& signals, const Dataset& data);


// a structure for passing information to minuit
struct Experiment {
  size_t nevents;
  std::vector<Signal>* signals;
  GPULL* ll;
};


// compute the negative log likelihood for the given parameters
// gpu-accelerated
void nll_gpu(int& ndim, double* gout, double& result, double* par, int flags);


// perform a fit or signal pdfs to the data, returning a TMinuit with final results
TMinuit* run_fit(const std::vector<Signal>& signals, const Dataset& data);

