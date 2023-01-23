#include "discrete-kernel.hpp"
#include <functional>
#include <cmath>

/* Discrete Kernel Constructor */
DiscreteKernel::DiscreteKernel(double aAlpha, double aBeta, function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&)> aKernelFunction)
{
    alpha = aAlpha;
    beta = aBeta;
    kernelFunction = aKernelFunction;
    kernelalphabeta = [](Vector3D& p, Vector3D& q, Vector3D& T)
    {
        Vector3D pmq = p - q;
        Vector3D projection = T ^ pmq;
    };
}
