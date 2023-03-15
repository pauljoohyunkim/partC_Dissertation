#include <iostream>
#include "solver.cuh"

/* (DO NOT USE DIRECTLY) Abstract Function: Addition of Tensor */
__global__ static void cuTensorAdd(CurveTensor& t1, CurveTensor& t2)
{
    int i = blockIdx.x;
    int j = blockIdx.y;
    int offset = i + (int) t1.N * j;

    if (offset < 3 * t1.N)
    {
        t1.dev_blocks[offset] = t1.dev_blocks[offset] + t2.dev_blocks[offset];
    }
}

CurveTensor::CurveTensor(unsigned int aN)
{
    N = aN;
    cudaMalloc((void**)&dev_blocks, 3 * N * sizeof(double));
    std::cout << "Tensor Constructed" << std::endl;
}

CurveTensor::~CurveTensor()
{
    cudaFree(dev_blocks);
    std::cout << "Tensor Destructed" << std::endl;
}

void tensorBlockLoad(CurveTensor& Gammabf, double* blocks, unsigned int N)
{
    cudaMemcpy(Gammabf.dev_blocks, blocks, 3 * N * sizeof(double), cudaMemcpyHostToDevice);
}
void tensorBlockFlush(CurveTensor& Gammabf, double* blocks, unsigned int N)
{
    cudaMemcpy(blocks, Gammabf.dev_blocks, 3 * N * sizeof(double), cudaMemcpyDeviceToHost);
}

void tensorAdd(CurveTensor& t1, CurveTensor& t2)
{
    dim3 grid(3, t1.N);
    cuTensorAdd<<<grid,1>>>(t1, t2);
}

