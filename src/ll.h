#ifndef __LL_H__
#define __LL_H__

#include <stdio.h>

// functor to encapsulate the cuda magic
class GPULL {
  public:
    GPULL(const float* lut, const size_t _nsignals, const size_t _nevents);
    virtual ~GPULL();
    void set_norms(const float* norms, const size_t _nsignals);
    double operator()(const double* _pars);

  private:
    size_t nsignals;
    size_t nevents;
    float* lut_device;
    float* norms_device;
    bool norms_device_alloc;
};

#endif  // __LL_H__

