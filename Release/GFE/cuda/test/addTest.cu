#include "../src/solver.cuh"

#define VECTORNUM 10

int main()
{
    CurveTensor tensor1 { VECTORNUM };
    CurveTensor tensor2 { VECTORNUM };
    double blocks1[3 * VECTORNUM] = { 1.0, 1.2, 1.3 };
    tensorBlockLoad(tensor1, blocks1, VECTORNUM);
    double blocks2[3 * VECTORNUM] = { 3.2, 4.2, -1.3 };
    tensorBlockLoad(tensor2, blocks2, VECTORNUM);

    //cuTensorAdd<<<grid,1>>> (tensor1, tensor2);
    tensorAdd(tensor1, tensor2);
    tensorBlockFlush(tensor1, blocks2, VECTORNUM);

    return 0;
    
}
