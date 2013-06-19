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
      delete[] Fit::norms;
      Fit::norms = nullptr;

      delete Fit::data;
      Fit::data = nullptr;

      delete this->minuit;
    }

    // run the fit, optionally specifying some different starting normalizations
    // and fixed/unfixed status
    TMinuit* operator()(Range<float> _e_range, Range<float> _r_range,
                        float live_time,
                        std::map<std::string, float>* _norms=nullptr,
                        std::map<std::string, bool>* _fix=nullptr,
                        bool run_minos=false);

  protected:
    // minuit fit function
    static void nll(int& ndim, double* gout, double& result, double* par,
                    int flags);

  private:
    TMinuit* minuit;
    static TH1* data;
    static float* norms;
    static std::vector<Signal> signals;
    static Range<float> e_range;
    static Range<float> r_range;
};

#endif  // __FIT_H__

