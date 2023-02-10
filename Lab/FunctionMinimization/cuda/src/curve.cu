#include "curve.hpp"
#include <vector>
#include <iostream>

/* Constructor */
FourierCurve::FourierCurve(std::vector<double>& axa, std::vector<double>& axb, std::vector<double>& aya, std::vector<double>& ayb, std::vector<double>& aza, std::vector<double>& azb, unsigned int aresolution)
{
    /* Pass initial coefficient list. */
    xa = axa;
    xb = axb;
    ya = aya;
    yb = ayb;
    za = aza;
    zb = azb;

    /* Set resolution of which will be used for computation. Higher resolution means more accurate result but lower performance. */
    resolution = aresolution;

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
    if (dev_trig_val_table_allocated)
    {
        cudaFree(dev_cos_table);
        cudaFree(dev_sin_table);
        std::cout << "dev_trig_val_table deallocated" << std::endl;
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
    if (!dev_trig_val_table_allocated)
    {
        /* Allocate VRAM */
        cudaMalloc((void**) &dev_cos_table, sizeof(double) * J * resolution);
        cudaMalloc((void**) &dev_sin_table, sizeof(double) * J * resolution);

        /* Precomputed Trig Value Generation
           This table takes the form of

                                 | cos 0t | cos 1t | cos 2t | ... | cos Jt 
         t = 0                   |
         t = 2pi/res * 1         |
         t = 2pi/res * 2         |
         ...                     |
         t = 2pi/res * (res-1)   |

           But flattened.
         */

        /* Filling the table row by row */
        std::vector<double> cos_row;
        std::vector<double> sin_row;
        cos_row.resize(J + 1);
        sin_row.resize(J + 1);
        for (unsigned int t = 0; t < resolution; t++)
        {
            for (unsigned int rowIndex = 0; rowIndex <= J; rowIndex++)
            {
                double theta = 2 * M_PI * t / resolution;
                cos_row[rowIndex] = cos(theta * rowIndex);
                sin_row[rowIndex] = sin(theta * rowIndex);
                cudaMemcpy(dev_cos_table + (J + 1) * t, &cos_row[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
                cudaMemcpy(dev_sin_table + (J + 1) * t, &sin_row[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
            }
        }

        dev_trig_val_table_allocated = true;
    }
}

__global__ void printCoefficientsPartiallyDEBUG(double* device_float_value)
{
    printf("%f\n", *device_float_value);
}
