#include "../src/solver.cuh"
#include <cstdio>

#define VECTORNUM 1

__global__ void kernel(double* dev_blocks)
{
    auto val = componentAccess(dev_blocks, 0, 2, VECTORNUM);
    printf("%f\n", val);
}

int main()
{
    CurveTensor tensor1 { VECTORNUM };
    double blocks1[3 * VECTORNUM] = { 1.0, 1.2, 1.3 };
    tensorBlockLoad(tensor1, blocks1);

    //cuTensorAdd<<<grid,1>>> (tensor1, tensor2);
    kernel<<<1,1>>> (tensor1.dev_blocks);

    return 0;
    
}
