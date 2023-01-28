#include "../src/solver.hpp"
#include <iostream>
#define N 100000

__global__ static void crossAll(double *dev_x1, double *dev_y1, double *dev_z1, double *dev_x2, double *dev_y2, double *dev_z2, double *dev_xv, double *dev_yv, double *dev_zv)
{
   int tid = blockIdx.x;
   if (tid < N)
   {
       cross(dev_x1[tid], dev_y1[tid], dev_z1[tid], dev_x2[tid], dev_y2[tid], dev_z2[tid], dev_xv[tid], dev_yv[tid], dev_zv[tid]);
   }
}

__global__ static void normAll(double *dev_x, double *dev_y, double *dev_z, double *dev_normv)
{
    int tid = blockIdx.x;
    if (tid < N)
    {
        dev_normv[tid] = l2norm3D(dev_x[tid], dev_y[tid], dev_z[tid]);
    }
}

int main()
{
    double x1[N];
    double y1[N];
    double z1[N];
    double x2[N];
    double y2[N];
    double z2[N];

    /* Filling in data */
    for (int i = 0; i < N; i++)
    {
        double val = (double) i;
        x1[i] = val * val;
        y1[i] = -val * val * 0.5;
        z1[i] = val;
        x2[i] = val * 0.7;
        y2[i] = -val;
        z2[i] = val * 0.3;
    }

    /* Results */
    double xv[N];
    double yv[N];
    double zv[N];
    double normsv[N];

    double *dev_x1, *dev_y1, *dev_z1, *dev_x2, *dev_y2, *dev_z2, *dev_xv, *dev_yv, *dev_zv, *dev_normv;

    /* Malloc on GPU */
    cudaMalloc((void**)&dev_x1, N * sizeof(double));
    cudaMalloc((void**)&dev_y1, N * sizeof(double));
    cudaMalloc((void**)&dev_z1, N * sizeof(double));
    cudaMalloc((void**)&dev_x2, N * sizeof(double));
    cudaMalloc((void**)&dev_y2, N * sizeof(double));
    cudaMalloc((void**)&dev_z2, N * sizeof(double));
    cudaMalloc((void**)&dev_xv, N * sizeof(double));
    cudaMalloc((void**)&dev_yv, N * sizeof(double));
    cudaMalloc((void**)&dev_zv, N * sizeof(double));
    cudaMalloc((void**)&dev_normv, N * sizeof(double));

    /* Copy the x, y, z array */
    cudaMemcpy(dev_x1, x1, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y1, y1, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_z1, z1, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_x2, x2, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y2, y2, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_z2, z2, N * sizeof(double), cudaMemcpyHostToDevice);

    /* Compute Cross Product for all */
    crossAll<<<N, 1>>>(dev_x1, dev_y1, dev_z1, dev_x2, dev_y2, dev_z2, dev_xv, dev_yv, dev_zv);

    /* Compute norm of the resultant cross products */
    normAll<<<N, 1>>>(dev_xv, dev_yv, dev_zv, dev_normv);

    /* Copy results to the result vector */
    cudaMemcpy(xv, dev_xv, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(yv, dev_yv, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(zv, dev_zv, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(normsv, dev_normv, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_x1);
    cudaFree(dev_y1);
    cudaFree(dev_z1);
    cudaFree(dev_x2);
    cudaFree(dev_y2);
    cudaFree(dev_z2);
    cudaFree(dev_xv);
    cudaFree(dev_yv);
    cudaFree(dev_zv);
    cudaFree(dev_normv);

    return 0;
}
