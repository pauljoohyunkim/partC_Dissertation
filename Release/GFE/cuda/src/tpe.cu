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

__device__ void ijk(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha, double beta)
{
    int i { p };
    int j { q };
    int k { r };
    Vector xkEdge = vectorFromTensor(dev_blocks, k+1, N) - vectorFromTensor(dev_blocks, k, N);
    double xkEdgeLen = xkEdge.norm();
    Vector xij = vectorFromTensor(dev_blocks, i, N) - vectorFromTensor(dev_blocks, j, N);
    double xijLen = xij.norm();

    xi = pow(xkEdgeLen * xijLen, 2) - pow(xkEdge * xij, 2);
    eta = pow(xijLen, beta) * pow(xkEdgeLen, alpha);
    dxi = -2 * xkEdge * pow(xijLen, 2) + 2 * (xkEdge * xij) * xij;
    deta = alpha * pow(xijLen, beta) * pow(xkEdgeLen,alpha-2) * (-xkEdge);
}

__device__ void km1jkm1(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha, double beta)
{
    int k { p + 1 };
    int j { q };
    Vector xkEdge = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, k-1, N);
    double xkEdgeLen = xkEdge.norm();
    Vector xkmj = vectorFromTensor(dev_blocks, k-1, N) - vectorFromTensor(dev_blocks, j, N);
    double xkmjLen = xkmj.norm();

    xi = pow(xkEdgeLen * xkmjLen, 2) - pow(xkEdge * xkmj, 2);
    eta = pow(xkmjLen, beta) * pow(xkEdgeLen, alpha);
    dxi = 2 * pow(xkmjLen, 2) * xkEdge - 2 * (xkEdge * xkmj) * xkmj;
    deta = alpha * pow(xkmjLen, beta) * pow(xkEdgeLen, alpha-2) * xkEdge;
}

__device__ void kjkm1(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha, double beta)
{
    int k { p };
    int j { q };
    Vector xkEdge = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, k-1, N);
    double xkEdgeLen = xkEdge.norm();
    Vector xkj = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, j, N);
    double xkjLen = xkj.norm();

    xi = pow(xkEdgeLen * xkjLen, 2) - pow(xkEdge * xkj, 2);
    eta = pow(xkjLen, beta) * pow(xkEdgeLen, alpha);
    dxi = 2 * pow(xkjLen, 2) * xkEdge + 2 * pow(xkEdgeLen, 2) * xkj
        - 2 * (xkEdge * xkj) * (xkEdge + xkj);
    deta = beta * pow(xkEdgeLen, alpha) * pow(xkjLen, beta-2) * xkj
        + alpha * pow(xkjLen, beta) * pow(xkEdgeLen, alpha-2) * xkEdge;
}
