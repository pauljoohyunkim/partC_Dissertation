#include <cstdio>
#include "solver.cuh"
#include "vector.cuh"
#include "tpe.cuh"

__device__ void dkalphabeta(double* dev_blocks, int p, int q, int r, int k, unsigned int N, Vector& res, double alpha, double beta)
{
    p = p % (int) N;
    q = q % (int) N;
    r = r % (int) N;
    k = k % (int) N;
    int km1 = ((k-1) % (int) N + N) % (int) N;
    double xi;
    double eta;
    Vector dxi;
    Vector deta;

    if (p == k && r == k)
    {
        kjk(dev_blocks, p, q, r, N, xi, eta, dxi, deta, alpha, beta);
        printf("kjk\n");
    }
    else if (r == k)
    {
        ijk(dev_blocks, p, q, r, N, xi, eta, dxi, deta, alpha, beta);
        printf("ijk\n");
    }
    else if (p == km1 && r == km1)
    {
        km1jkm1(dev_blocks, p, q, r, N, xi, eta, dxi, deta, alpha, beta);
        printf("km1jkm1\n");
    }
    else if (p == k && r == km1)
    {
        kjkm1(dev_blocks, p, q, r, N, xi, eta, dxi, deta, alpha, beta);
        printf("kjkm1\n");
    }
    else if (q == k)
    {
        ikj(dev_blocks, p, q, r, N, xi, eta, dxi, deta, alpha, beta);
        printf("ikj\n");
    }
    else
    {
        printf("(p,q,r) tuple not defined\n");
    }

    res = (alpha / 2 * pow(xi, alpha/2-1) * dxi * eta - pow(xi, alpha/2) * deta) / pow(eta,2);
}

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

__device__ void ikj(double* dev_blocks, int p, int q, int r, unsigned int N, double& xi, double&eta, Vector& dxi, Vector& deta, double alpha, double beta)
{
    int i { p };
    int j { r };
    int k { q };
    Vector xjEdge = vectorFromTensor(dev_blocks, j+1, N) - vectorFromTensor(dev_blocks, j, N);
    double xjEdgeLen = xjEdge.norm();
    Vector xki = vectorFromTensor(dev_blocks, k, N) - vectorFromTensor(dev_blocks, i, N);
    double xkiLen = xki.norm();

    xi = pow(xjEdgeLen * xkiLen, 2) - pow(xjEdge * xki, 2);
    eta = pow(xkiLen, beta) * pow(xjEdgeLen, alpha);
    dxi = 2 * pow(xjEdgeLen, 2) * xki - 2 * (xjEdge * xki) * xjEdge;
    deta = beta * pow(xjEdgeLen, alpha) * pow(xkiLen, beta-2) * xki;
}
