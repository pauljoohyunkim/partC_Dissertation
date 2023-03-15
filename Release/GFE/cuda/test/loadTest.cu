#include "../src/solver.cuh"

#define VECTORNUM 10

/*
__global__ void testKernel(double* blockAddr1, double* blockAddr2, unsigned int N)
{
    cuCurveTensor tensor1 { N };
}
*/

int main()
{
    CurveTensor tensor1 { VECTORNUM };
    double blocks1[3 * VECTORNUM] = { 1.0, 1.2, 1.3 };
    tensorBlockLoad(tensor1, blocks1, VECTORNUM);
    double blocks2[3 * VECTORNUM] = { 0 };
    tensorBlockFlush(tensor1, blocks2, VECTORNUM);

    return 0;
    
}
