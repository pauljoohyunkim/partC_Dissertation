#include <vector>
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
    /* Constructor Test 1 */
    CurveTensor tensor1 { VECTORNUM };
    double blocks1[3 * VECTORNUM] = { 1.0, 1.2, 1.3 };
    tensorBlockLoad(tensor1, blocks1);
    double blocks2[3 * VECTORNUM] = { 0 };
    tensorBlockFlush(tensor1, blocks2);

    /* Constructor Test 2 */
    std::vector<double> x { 1.0, 2.0, 3.0, 4.0, 5.0 };
    std::vector<double> y { -2.0, 2.3, -4.2, -30, 2 };
    std::vector<double> z { 3.0, 4.3, -1.2, 0, 2.2 };
    CurveTensor tensor3 { x, y, z };
    double blocks3[3 * VECTORNUM] = { 0 };
    tensorBlockFlush(tensor3, blocks3);

    return 0;
    
}
