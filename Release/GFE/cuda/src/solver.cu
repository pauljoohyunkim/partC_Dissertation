#include <iostream>
#include <vector>
#include "solver.cuh"

/* (DO NOT USE DIRECTLY) Abstract Function: Addition of Tensor */
__global__ static void cuTensorAdd(double* dev_blocks_1, double* dev_blocks_2, int N, double scalar=1)
{
    int i = blockIdx.x;
    int j = blockIdx.y;
    int offset = N * i + j;

    if (offset < 3 * N)
    {
        dev_blocks_1[offset] += scalar * dev_blocks_2[offset];
    }
}

/* (DO NOT USE DIRECTLY) Abstract Function: Subtraction of Tensor */
__global__ static void cuTensorSubtract(double* dev_blocks_1, double* dev_blocks_2, int N, double scalar=1)
{
    int i = blockIdx.x;
    int j = blockIdx.y;
    int offset = N * i + j;

    if (offset < 3 * N)
    {
        dev_blocks_1[offset] -= scalar * dev_blocks_2[offset];
    }
}

CurveTensor::CurveTensor(unsigned int aN)
{
    N = aN;
    cudaMalloc((void**)&dev_blocks, 3 * N * sizeof(double));
    std::cout << "Tensor Constructed" << std::endl;
}

CurveTensor::CurveTensor(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{
    N = x.size();
    cudaMalloc((void**)&dev_blocks, 3 * N * sizeof(double));
    cudaMemcpy(dev_blocks, &x[0], N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_blocks + N, &y[0], N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_blocks + 2 * N, &z[0], N * sizeof(double), cudaMemcpyHostToDevice);
}

CurveTensor::~CurveTensor()
{
    cudaFree(dev_blocks);
    std::cout << "Tensor Destructed" << std::endl;
}

__device__ double& componentAccess(double* dev_blocks, int i, unsigned int j, unsigned int N)
{
    i = (i % (int) N + N) % (int) N;
    return dev_blocks[N * (unsigned int) j + i];
}

void tensorBlockLoad(CurveTensor& Gammabf, double* blocks)
{
    unsigned int N { Gammabf.N };
    cudaMemcpy(Gammabf.dev_blocks, blocks, 3 * N * sizeof(double), cudaMemcpyHostToDevice);
}
void tensorBlockFlush(CurveTensor& Gammabf, double* blocks)
{
    unsigned int N { Gammabf.N };
    cudaMemcpy(blocks, Gammabf.dev_blocks, 3 * N * sizeof(double), cudaMemcpyDeviceToHost);
}
void tensorBlockFlush(CurveTensor& Gammabf, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{
    unsigned int N { Gammabf.N };
    //cudaMemcpy(blocks, Gammabf.dev_blocks, 3 * N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(x[0]), Gammabf.dev_blocks, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(y[0]), Gammabf.dev_blocks + N, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(z[0]), Gammabf.dev_blocks + 2 * N, N * sizeof(double), cudaMemcpyDeviceToHost);
}

void tensorAdd(CurveTensor& t1, CurveTensor& t2, double scalar)
{
    dim3 grid(3, t1.N);
    cuTensorAdd<<<grid,1>>>(t1.dev_blocks, t2.dev_blocks, (int) t1.N, scalar);
}

void tensorSubtract(CurveTensor& t1, CurveTensor& t2, double scalar)
{
    dim3 grid(3, t1.N);
    cuTensorSubtract<<<grid,1>>>(t1.dev_blocks, t2.dev_blocks, (int) t1.N, scalar);
}

