#include <iostream>
#include "solver.cuh"

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
