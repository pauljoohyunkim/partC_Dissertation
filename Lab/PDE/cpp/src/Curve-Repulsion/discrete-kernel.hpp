#ifndef DISCRETE_KERNEL_HPP
#define DISCRETE_KERNEL_HPP

#include "../math-objects.hpp"
#include <functional>

class DiscreteKernel
{
    public:
        /* Discrete Kernel Constructor */
        /* For kernel function k_{ij}, arguments are x_i, x_{i+1}, x_j, x_{j+1}, T_{i} */
        DiscreteKernel(double aAlpha, double aBeta, function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&)> aKernelFunction);

        function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&)> kernelFunction;
    private:
        double alpha {};
        double beta {};
        function<double(Vector3D&, Vector3D&, Vector3D&)> kernelalphabeta;
};


#endif  // discrete-kernel.hpp
