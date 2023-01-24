#include "discrete-kernel.hpp"
#include <functional>
#include <cmath>
#include <cstdlib>

/* Discrete Kernel Constructor */
DiscreteKernel::DiscreteKernel(double aAlpha, double aBeta, std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, DiscreteKernel&)> aKernelFunction)
{
    alpha = aAlpha;
    beta = aBeta;
    kernelFunction = aKernelFunction;
    kernelalphabeta = [this](Vector3D& p, Vector3D& q, Vector3D& T)
    {
        Vector3D pmq = p - q;
        Vector3D projection = T ^ pmq;
        double numerator = pow(l2norm(projection), alpha);
        double denominator = pow(l2norm(pmq), beta);
        return numerator / denominator;
    };
}

/* Discretized Energy */
double DiscreteKernel::energy(Curve C)
{
    double e { 0 };
    auto J = (int) C.getNPoints();

    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J ; j++)
        {
            if (abs(i - j) > 1 && abs(i - j + J) > 1 && abs(i - j - J) > 1)
            {
                /* Components of E
                 * p = x_i
                 * q = x_j
                 * pI = x_{i+1} - x_i
                 * lI = |pI|
                 * TI = pI / lI ("Tangent")
                 * */
                auto p = C[i];
                auto q = C[j];
                auto pI = C[i+1] - p;
                auto lI = l2norm(pI);
                auto qJ = C[j+1] - q;
                auto lJ = l2norm(qJ);
                auto TI = pI * (1 / lI);
                e += kernelFunction(p, C[i + 1], q, C[j + 1], TI, *this) * lI * lJ;
            }
        }
    }

    return e;

}
