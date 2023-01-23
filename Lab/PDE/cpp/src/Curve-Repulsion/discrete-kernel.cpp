#include "discrete-kernel.hpp"
#include <functional>
#include <cmath>

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
