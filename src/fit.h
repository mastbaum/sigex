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
class Dataset;


// user interface to the MINUIT ML fit
//
// WARNING: HEREIN LIE GLOBAL MEMBER VARIABLES. MINUIT USES GLOBAL STATE.
//          YOU CANNOT DO MORE THAN ONE FIT AT A TIME.
class Fit {
  public:
    Fit(const std::vector<Signal>& signals, const Dataset& data);

    virtual ~Fit() {
      delete Fit::ll;
      Fit::ll = nullptr;

      delete this->lut;
      delete this->minuit;
    }

    // run the fit, optionally specifying some different starting normalizations
    // and fixed/unfixed status
    TMinuit* operator()(Range<float> e_range, Range<float> r_range,
                        std::map<std::string, float>* _norms=nullptr,
                        std::map<std::string, bool>* _fix=nullptr,
                        bool run_minos=false);

  protected:
    // construct a P(x) lookup table
    static float* build_lut(const std::vector<Signal>& signals, const Dataset& data);

    // minuit fit function
    static void nll(int& ndim, double* gout, double& result, double* par,
                    int flags);

  private:
    float* lut;  // P(x) lookup table
    TMinuit* minuit;
    static size_t nevents;  // total event counts in dataset (rows in lut)
    static GPULL* ll;  // GPU LL interface
    static std::vector<Signal> signals;
};

#endif  // __FIT_H__

