#include <cuda.h>
#include "ll.h"

// macro to print error and abort on cuda errors, c/o stan
#define CUDA_CHECK_ERROR(call) do { \
  cudaError err = call; \
  if (cudaSuccess != err) { \
    fprintf(stderr, "Cuda error in file '%s' in line %i: %s.\n", \
            __FILE__, __LINE__, cudaGetErrorString(err)); \
    exit(EXIT_FAILURE); \
  } \
} while (0)


// compute the log likelihood given signal rates and normalizations
__global__ void ll(const float* lut, const float* n, const double* pars, const size_t ne, const size_t ns, double* sums) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  for (int i=idx; i<(int)ne; i+=gridDim.x*blockDim.x) {
    double s = 0;
    for (size_t j=0; j<ns; j++) {
      s += pars[j] * n[j] * lut[i*ns+j];
    }
    sums[idx] += log(s);
  }
}


GPULL::GPULL(const float* lut, const size_t _nsignals, const size_t _nevents)
  : nsignals(_nsignals), nevents(_nevents), norms_device_alloc(false) {
  CUDA_CHECK_ERROR(cudaMalloc(&this->lut_device, this->nevents * this->nsignals * sizeof(float)));
  CUDA_CHECK_ERROR(cudaMemcpy(this->lut_device, lut, this->nevents * this->nsignals * sizeof(float), cudaMemcpyHostToDevice));
}


GPULL::~GPULL() {
  CUDA_CHECK_ERROR(cudaFree(this->norms_device));
  CUDA_CHECK_ERROR(cudaFree(this->lut_device));
}


void GPULL::set_norms(const float* norms, const size_t _nsignals) {
  if (!this->norms_device_alloc) {
    CUDA_CHECK_ERROR(cudaMalloc(&this->norms_device, this->nsignals * sizeof(float)));
    this->norms_device_alloc = true;
  }
  CUDA_CHECK_ERROR(cudaMemcpy(this->norms_device, norms, this->nsignals * sizeof(float), cudaMemcpyHostToDevice));
}


double GPULL::operator()(const double* pars) {
  size_t blocksize = 256;
  size_t nblocks = 16;
  size_t nthreads = nblocks * blocksize;

  double* pars_device;
  CUDA_CHECK_ERROR(cudaMalloc(&pars_device, this->nsignals * sizeof(double)));
  CUDA_CHECK_ERROR(cudaMemcpy(pars_device, pars, this->nsignals * sizeof(double), cudaMemcpyHostToDevice));

  double* sums_device;
  CUDA_CHECK_ERROR(cudaMalloc(&sums_device, nthreads * sizeof(double)));
  CUDA_CHECK_ERROR(cudaMemset(sums_device, 0, nthreads * sizeof(double)));

  ll<<<nblocks, blocksize>>>(this->lut_device, this->norms_device, pars_device, this->nevents, this->nsignals, sums_device);
  CUDA_CHECK_ERROR(cudaThreadSynchronize());

  double* sums = new double[blocksize * nblocks];
  CUDA_CHECK_ERROR(cudaMemcpy(sums, sums_device, nthreads * sizeof(double), cudaMemcpyDeviceToHost));

  CUDA_CHECK_ERROR(cudaFree(pars_device));
  CUDA_CHECK_ERROR(cudaFree(sums_device));

  double sum = 0;
  for (size_t i=0; i<nthreads; i++) {
    sum += sums[i];
  }

  delete[] sums;

  return sum;
}

