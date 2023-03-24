#include <cstdio>
#include "solver.cuh"
#include "vector.cuh"
#include "tpe.cuh"

template <class T>
ScratchPad<T>::ScratchPad(unsigned int am, unsigned int aN)
{
    m = am;
    N = aN;
    cudaMalloc((void**)(&scratchpads), m * N * sizeof(T));
}

template <class T>
ScratchPad<T>::~ScratchPad()
{
    cudaFree(scratchpads);
}




/* Static Function Declaration */
//__global__ static void cuDEnergy(double* dev_curve_tensor_blocks, double* dev_differential_blocks, int** derivative_index_scratch, unsigned int N, double alpha=3, double beta=6);




__device__ double kernelalphabeta(Vector& p, Vector& q, Vector& T, double alpha, double beta)
{
    Vector pmq = p - q;
    auto numerator = pow(norm(T ^ pmq), alpha);
    auto denominator = pow(pmq.norm(), beta);

    return numerator / denominator;
}

__device__ double kij(double* dev_blocks, int i, int j, unsigned int N, double alpha, double beta)
{
    Vector TI = vectorFromTensor(dev_blocks, i+1, N) - vectorFromTensor(dev_blocks, i, N);
    TI = TI / TI.norm();
    
    double res = 0;
    Vector xi = vectorFromTensor(dev_blocks, i, N);
    Vector xj = vectorFromTensor(dev_blocks, j, N);
    Vector xip = vectorFromTensor(dev_blocks, i+1, N);
    Vector xjp = vectorFromTensor(dev_blocks, j+1, N);
    res += kernelalphabeta(xi, xj, TI, alpha, beta);
    res += kernelalphabeta(xi, xjp, TI, alpha, beta);
    res += kernelalphabeta(xip, xj, TI, alpha, beta);
    res += kernelalphabeta(xip, xjp, TI, alpha, beta);

    return res;
}

__device__ void dkij(double* dev_blocks, int i, int j, int k, unsigned int N, Vector& res, double alpha, double beta)
{
    res = Vector();
    Vector temp;

    dkalphabeta(dev_blocks, i, j, i, k, N, temp, alpha, beta);
    res = res + temp;
    dkalphabeta(dev_blocks, i, j+1, i, k, N, temp, alpha, beta);
    res = res + temp;
    dkalphabeta(dev_blocks, i+1, j, i, k, N, temp, alpha, beta);
    res = res + temp;
    dkalphabeta(dev_blocks, i+1, j+1, i, k, N, temp, alpha, beta);
    res = res + temp;
    res = res / 4;
}

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
        res = Vector(0, 0, 0);
        return;
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


__device__ void dProductOfLengths(double* dev_blocks, int p, int q, int k, unsigned int N, Vector& res)
{
    int i { p % (int) N };
    int j { q % (int) N };
    k = k % (int) N;
    int ip1 { (i+1) % (int) N };
    int jp1 { (j+1) % (int) N };

    Vector xiEdge = vectorFromTensor(dev_blocks, i+1, N) - vectorFromTensor(dev_blocks, i, N);
    double xiEdgeLen = xiEdge.norm();
    Vector xjEdge = vectorFromTensor(dev_blocks, j+1, N) - vectorFromTensor(dev_blocks, j, N);
    double xjEdgeLen = xjEdge.norm();

    if (i == k)
    {
        res = xjEdgeLen * (-xiEdge) / xiEdgeLen;
    }
    else if (ip1 == k)
    {
        res = xjEdgeLen * xiEdge / xiEdgeLen;
    }
    else if (j == k)
    {
        res = xiEdgeLen * (-xjEdge) / xjEdgeLen;
    }
    else if (jp1 == k)
    {
        res = xiEdgeLen * xjEdge / xjEdgeLen;
    }
    else
    {
        res = Vector(0, 0, 0);
    }
}

/* 4(N-1) Pairs of Form:
   * *
   * *
   * *
 */
__device__ void fillDerivativeIndex(int* dev_derivative_indices, int k, unsigned int J)
{
    int N { (int) J };
    int stackpt = 0;
    for (int i = 0; i < N; i++)
    {
        if (abs(k-i) > 1 && abs(k-i+N) > 1 && abs(k-i-N) > 1)
        {
            dev_derivative_indices[2 * stackpt] = k;
            dev_derivative_indices[2 * stackpt + 1] = i;
            stackpt++;
        }
    }
    for (int i = 0; i < N; i++)
    {
        if (abs(k-i) > 1 && abs(k-i+N) > 1 && abs(k-i-N) > 1)
        {
            dev_derivative_indices[2 * stackpt] = i;
            dev_derivative_indices[2 * stackpt + 1] = k;
            stackpt++;
        }
    }
    for (int i = 0; i < N; i++)
    {
        if (abs(k-1-i) > 1 && abs(k-1-i+N) > 1 && abs(k-1-i-N) > 1)
        {
            dev_derivative_indices[2 * stackpt] = (k+N-1)%N;
            dev_derivative_indices[2 * stackpt + 1] = i;
            stackpt++;
        }
    }
    for (int i = 0; i < N; i++)
    {
        if (abs(k-1-i) > 1 && abs(k-1-i+N) > 1 && abs(k-1-i-N) > 1)
        {
            dev_derivative_indices[2 * stackpt] = i;
            dev_derivative_indices[2 * stackpt + 1] = (k+N-1)%N;
            stackpt++;
        }
    }
}

/* <<<N, 1>>> 

 Note that derivative_index_scratch should be called from ScratchPad<int> { N, 8 * (N - 3) }*/
//__global__ static void cuDEnergy(double* dev_curve_tensor_blocks, double* dev_differential_blocks, int** derivative_index_scratch, unsigned int N, double alpha, double beta)
__global__ void cuDEnergy(double* dev_curve_tensor_blocks, double* dev_differential_blocks, int* derivative_index_scratch, unsigned int N, double alpha, double beta)
//__global__ static void cuDEnergy(double* dev_curve_tensor_blocks, double* dev_differential_blocks, int** derivative_index, unsigned int N, double alpha, double beta)
{
    int k = blockIdx.x;

    /* This vector will be the column of the tensor */
    Vector res(0, 0, 0);

    /* Generate the indices relevant to derivative for each k */
    int* dev_derivative_indices = derivative_index_scratch + N * k;
    //int dev_derivative_indices [8 * (N - 3)];
    fillDerivativeIndex(dev_derivative_indices, k, N);

    /* Loop over each pair of indices */
    for (int index = 0; index < 4 * ((int)N-3); index++)
    {
        int i = dev_derivative_indices[2 * index];
        int j = dev_derivative_indices[2 * index + 1];

        Vector xiEdge = vectorFromTensor(dev_curve_tensor_blocks, i+1, N) - vectorFromTensor(dev_curve_tensor_blocks, i, N);
        double xiEdgeLen = xiEdge.norm();
        Vector xjEdge = vectorFromTensor(dev_curve_tensor_blocks, j+1, N) - vectorFromTensor(dev_curve_tensor_blocks, j, N);
        double xjEdgeLen = xjEdge.norm();

        /* Each summand is from the product rule of k_ij and product of edge lengths */
        Vector summand1;
        Vector summand2;
        dkij(dev_curve_tensor_blocks, i, j, k, N, summand1, alpha, beta);
        summand1 = summand1 * xiEdgeLen * xjEdgeLen;

        dProductOfLengths(dev_curve_tensor_blocks, i, j, k, N, summand2);
        summand2 = kij(dev_curve_tensor_blocks, i, j, N, alpha, beta) * summand2;
        
        res = res + summand1;
        res = res + summand2;
    }

    componentAccess(dev_differential_blocks, k, 0, N) = res.x;
    componentAccess(dev_differential_blocks, k, 1, N) = res.y;
    componentAccess(dev_differential_blocks, k, 2, N) = res.z;
}




/* Explicit Instantiation */
template class ScratchPad<int>;
