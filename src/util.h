#ifndef __UTIL_H__
#define __UTIL_H__

// Various utilities
#include <vector>
#include <string>
#include <assert.h>
#include <string>
#include <TTree.h>
#include <TChain.h>

class TH1;
class TH1D;
class TH2F;
class TNtuple;

// read a list of files out of one or more text files
std::vector<std::string> read_file_list(std::vector<std::string> l);


// make a tchain from the filenames, checking each one first to make sure it's
// not a zombie.
TChain* make_tchain(std::vector<std::string> filenames, std::string name="T");


// find the median value in a vector
double median(std::vector<double> v);


// container for a range of values
template <class T>
struct Range {
  T min;
  T max;
};


// container for signal metadata and pdfs
struct Signal {
  std::string name;
  std::string title;  // histogram title in ROOT-LaTeX format
  float rate;  // expected events per year
  float constraint;  // fractional uncertainty
  bool fixed;  // float normalization?
  int ndim;  // number of dimensions in histogram
  TH1* histogram;  // TH1F or TH2F -- see ndim
};


// types (/dimensionalities) of fit
enum class FitMode {
  ENERGY,
  ENERGY_RADIUS
};


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


// manages the configuration of the fit
// loads and parses a json config file
class FitConfig {
  public:
    FitConfig(std::string const filename);
    virtual ~FitConfig() {}

    // pretty-print the fit parameters
    void print() const;

    FitMode mode;  // fit type (energy, energy+radius, ...)
    float confidence;  // confidence level for results (e.g. 0.9)
    float live_time;  // experiment live time in years
    float efficiency;  // overall efficiency correction
    unsigned mc_trials;  // number of experiments to generate
    unsigned fake_experiments;  // number of fake experiments per experiment
    std::string output_file;  // base filename for output
    std::string signal_name;  // name of the signal that is the signal
    Range<float> e_range;  // range of energies to include in fit
    Range<float> r_range;  // range of radii to include in fit
    std::vector<Signal> signals;  // signal histogram and metadata

    // project a TH2F down to an energy-only TH1F, optionally cutting on radius
    static TH1D* project1d(TH2F* const hist2d, Range<float>* const r_range=nullptr);

  protected:
    // load an energy/radius ROOT TH2F from a file
    static TH2F* load_histogram(std::string const filename, std::string const objname);
};

#endif  // __UTIL_H__

