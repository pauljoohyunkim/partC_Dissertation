#ifndef __TPE_CUH__
#define __TPE_CUH__

#include "vector.cuh"

__device__ void kjk(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void ijk(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void km1jkm1(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);
__device__ void kjkm1(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha=3, double beta=6);

#endif  // __TPE_CUH__
