#ifndef DISCRETE_KERNEL_HPP
#define DISCRETE_KERNEL_HPP

#include "../math-objects.hpp"
#include <functional>

class DiscreteKernel
{
    public:
        /* Discrete Kernel Constructor */
        /* For kernel function k_{ij}, arguments are x_i, x_{i+1}, x_j, x_{j+1}, T_{i}, and the Discrete Kernel object itself */
        DiscreteKernel(double aAlpha, double aBeta, std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, DiscreteKernel&)> aKernelFunction);

        std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, DiscreteKernel&)> kernelFunction;
        std::function<double(Vector3D&, Vector3D&, Vector3D&)> kernelalphabeta;
    private:
        double alpha {};
        double beta {};
};


#endif  // discrete-kernel.hpp
