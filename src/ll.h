#include <stdio.h>

// functor to encapsulate the cuda magic
class GPULL {
  public:
    GPULL(const float* norms, const float* lut, const size_t _nsignals, const size_t _nevents);
    virtual ~GPULL();
    double operator()(const double* _pars);

  private:
    size_t nsignals;
    size_t nevents;
    float* lut_device;
    float* norms_device;
};

