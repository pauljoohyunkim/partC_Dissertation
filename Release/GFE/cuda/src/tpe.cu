#include "solver.cuh"
#include "vector.cuh"
#include "tpe.cuh"

__device__ void kjk(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha, double beta)
{
    int k { p };
    int j { q };
    Vector xkEdge = vectorFromTensor(dev_blocks, k+1, N) - vectorFromTensor(dev_blocks, k, N);
    double xkEdgeLen = xkEdge.norm();
    Vector xkj = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, j, N);
    double xkjLen = xkj.norm();

    xi = pow(xkEdgeLen * xkjLen, 2) - pow(xkEdge * xkj, 2);
    eta = pow(xkjLen, beta) * pow(xkEdgeLen, alpha);
    dxi = -2 * xkEdge * pow(xkjLen, 2) + 2 * pow(xkEdgeLen, 2) * xkj 
        - 2 * (xkEdge * xkj) * (xkEdge - xkj);
    deta = beta * pow(xkjLen, beta-2) * pow(xkEdgeLen, alpha) * xkj
        + alpha * pow(xkjLen, beta) * pow(xkEdgeLen, alpha-2) * (-xkEdge);
}
