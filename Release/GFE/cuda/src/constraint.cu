#include "vector.cuh"
#include "solver.cuh"
#include "tpe.cuh"

__global__ static void cuInitialEdgeLengthConstraintInit(double* dev_blocks, double* initEdgeLengths, unsigned int N)
{
    int i = blockIdx.x;

    /* Compute L_i for each i */
    if (i < (int) N)
    {
        Vector xiEdge = vectorFromTensor(dev_blocks, i+1, N) - vectorFromTensor(dev_blocks, i, N);
        initEdgeLengths[i] = xiEdge.norm();
    }

}
__global__ static void cuInitialEdgeLengthConstraint(double* dev_blocks, double* dev_d_constraint, double* initEdgeLengths, unsigned int N)
{
    int k = blockIdx.x;
    int km1 = ((k-1) % (int) N + N) % (int) N;

    /* Compute dC/dx_i for each i */
    if (k < (int) N)
    {
        xkEdge = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, k+1, N);
        xkEdgeNorm = xkEdge.norm();
        xkmEdge = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, km1, N);
        xkmEdgeNorm = xkmEdge.norm();

        Vector summand1;
        summand1 = 2 * (xkEdgeNorm - initEdgeLengths[k]) / xkEdgeNorm * xkEdge;
        Vector summand2;
        summand2 = 2 * (xkmEdgeNorm - initEdgeLengths[km1]) / xkmEdgeNorm * xkmEdge;
        Vector res = summand1 + summand2;

        componentAccess(dev_d_constraint, k, 0, N) = res.x;
        componentAccess(dev_d_constraint, k, 1, N) = res.y;
        componentAccess(dev_d_constraint, k, 2, N) = res.z;
    }
}

void initialEdgeLengthConstraint(CurveTensor& curve, CurveTensor& dConstraint, ScratchPad<double>& initEdgeLengths)
{
    /* Fill the initial edge lengths onto scratch pad of size 3*N */
    static bool init { true };
    unsigned int N = curve.N;
    if (init)
    {
        cuInitialEdgeLengthConstraintInit<<<N, 1>>>(curve.dev_blocks, dConstraint.dev_blocks, N);
        cudaDeviceSynchronize();
        init = false;
    }
    cuInitialEdgeLengthConstraint<<<N, 1>>>(curve.dev_blocks, dConstraint.dev_blocks, N);
    cudaDeviceSynchronize();
}
