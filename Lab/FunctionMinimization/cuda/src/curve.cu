#include "curve.hpp"
#include <vector>
#include <iostream>

/* Constructor */
FourierCurve::FourierCurve(std::vector<double>& axa, std::vector<double>& axb, std::vector<double>& aya, std::vector<double>& ayb, std::vector<double>& aza, std::vector<double>& azb)
{
    /* Pass initial coefficient list. */
    xa = axa;
    xb = axb;
    ya = aya;
    yb = ayb;
    za = aza;
    zb = azb;

    /* Up to j^th order term */
    J = axa.size() - 1;
}

/* Destructor */
FourierCurve::~FourierCurve()
{
    if (dev_coefficient_allocated)
    {
        cudaFree(dev_coefficients);
        std::cout << "dev_coefficients deallocated" << std::endl;
    }
    if (dev_scratch_pad_allocated)
    {
        cudaFree(dev_scratch_pad);
        std::cout << "dev_scratch_pad deallocated" << std::endl;
    }
}

/* Load onto GPU */
void FourierCurve::cudafy()
{
    if (!dev_coefficient_allocated)
    {
        /* Allocate VRAM */
        cudaMalloc((void**) &dev_coefficients, sizeof(double) * 6 * (J + 1));

        /* Copy all the coefficients, concatenated. */
        cudaMemcpy(dev_coefficients, &xa[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coefficients + (J + 1), &xb[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coefficients + 2 * (J + 1), &ya[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coefficients + 3 * (J + 1), &yb[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coefficients + 4 * (J + 1), &za[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coefficients + 5 * (J + 1), &zb[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
        dev_coefficient_allocated = true;
    }
}

__global__ void printCoefficientsPartiallyDEBUG(double* devcoeffs, int index)
{
    printf("%f\n", devcoeffs[index]);
}
