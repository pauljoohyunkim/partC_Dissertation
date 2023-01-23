#include "discrete-kernel.hpp"
#include <functional>

/* Discrete Kernel Constructor */
DiscreteKernel::DiscreteKernel(double aAlpha, double aBeta, function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&)> aKernelFunction)
{
    alpha = aAlpha;
    beta = aBeta;
    kernelFunction = aKernelFunction;
}
