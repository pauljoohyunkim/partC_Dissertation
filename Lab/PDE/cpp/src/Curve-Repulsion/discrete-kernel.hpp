#ifndef DISCRETE_KERNEL_HPP
#define DISCRETE_KERNEL_HPP

#include "../math-objects.hpp"
#include "../geometric-objects.hpp"
#include <functional>

class DiscreteKernel
{
    public:
        /* Discrete Kernel Constructor */
        /* For kernel function k_{ij}, arguments are x_i, x_{i+1}, x_j, x_{j+1}, T_{i}, and the Discrete Kernel object itself */
        DiscreteKernel(double aAlpha, double aBeta, std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, DiscreteKernel&)> aKernelFunction);

        /* Discretized Energy */
        double energy(Curve C);
        double energyDifferential(Curve C, Vector3D v, int index);



        /* k_{i,j} */
        std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, DiscreteKernel&)> kernelFunction;
        /* k_\beta^\alpha (Actual) */
        std::function<double(Vector3D&, Vector3D&, Vector3D&)> kernelalphabeta;

    private:
        double alpha { 2 };
        double beta { 4 };
};


#endif  // discrete-kernel.hpp
