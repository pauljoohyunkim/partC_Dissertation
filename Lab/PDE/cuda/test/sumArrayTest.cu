#include "../src/solver.hpp"
#include <vector>
#include <iostream>

#define LENGTH 4

__global__ void kernel(double* dev_x)
{
    double sum = sumArray(dev_x, LENGTH);

    dev_x[0] = sum;
}

int main()
{
    double arrOfVal[] = { 1, 2, 3, 4 };
    double* dev;
    double sum;

    /* Alloc */
    cudaMalloc((void**)&dev, sizeof(double) * LENGTH);
    std::cout << "cudaMalloc-ed" << std::endl;

    cudaMemcpy(dev, arrOfVal, sizeof(double) * LENGTH, cudaMemcpyHostToDevice);
    kernel<<<1, 1>>>(dev);
    cudaMemcpy(&sum, dev, sizeof(double), cudaMemcpyDeviceToHost);

    std::cout << sum << std::endl;
    
    cudaFree(dev);
}
