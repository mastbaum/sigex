#ifndef __UTIL_H__
#define __UTIL_H__

// Various utilities
#include <vector>
#include <string>
#include <assert.h>
#include <string>

class TChain;
class TH1;
class TH1D;
class TH2F;

// read a list of files out of one or more text files
std::vector<std::string> read_file_list(std::vector<std::string> l);


// make a tchain from the filenames, checking each one first to make sure it's
// not a zombie.
TChain* make_tchain(std::vector<std::string> filenames, std::string name="T");


// find the median value in a vector
double median(std::vector<double> &v);


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
    unsigned mc_trials;  // number of fake experiments to generate
    std::string signal_name;  // name of the signal that is the signal
    Range<float> e_range;  // range of energies to include in fit
    Range<float> r_range;  // range of radii to include in fit
    std::vector<Signal> signals;  // signal histogram and metadata

  protected:
    // load an energy/radius ROOT TH2F from a file
    static TH2F* load_histogram(std::string const filename, std::string const objname);

    // project a TH2F down to an energy-only TH1F, optionally cutting on radius
    static TH1D* project1d(TH2F* const hist2d, Range<float>* const r_range=nullptr);
};

#endif  // __UTIL_H__

