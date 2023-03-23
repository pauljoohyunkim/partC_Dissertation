#ifndef __TPE_CUH__
#define __TPE_CUH__

#include "vector.cuh"

template <class T>
struct ScratchPad
{
    /* m arrays of length N */
    ScratchPad(unsigned int am, unsigned int aN);
    ~ScratchPad();

    T** scratchpads;
    unsigned int m;
    unsigned int N;

};

__device__ double kernelalphabeta(Vector& p, Vector& q, Vector& T, double alpha=3, double beta=6);
__device__ double kij(double* dev_blocks, int i, int j, unsigned int N, double alpha=3, double beta=6);

__device__ void dkij(double* dev_blocks, int i, int j, int k, unsigned int N, Vector& res, double alpha=3, double beta=6);
__device__ void dkalphabeta(double* dev_blocks, int p, int q, int r, int k, unsigned int N, Vector& res, double alpha=3, double beta=6);
__device__ void kjk(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void ijk(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void km1jkm1(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void kjkm1(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void ikj(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);

__device__ void dProductOfLengths(double* dev_blocks, int p, int q, int k, unsigned int N, Vector& res);

__device__ void fillDerivativeIndex(int* dev_derivative_indices, int k, unsigned int J);

#endif  // __TPE_CUH__
