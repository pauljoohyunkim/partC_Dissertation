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

    /* Set host memory for curve points */
    x.resize(resolution, 0);
    y.resize(resolution, 0);
    z.resize(resolution, 0);

    /* differential */
    coeff_differential.resize(6 * (J + 1), 0);

    /* energyMatrix (For DEBUG purposes) */
    energyMatrix.resize(resolution * resolution, 0);

    std::cout << "Fourier curve constructed with" << std::endl;
    std::cout << "J: " << J << std::endl;
    std::cout << "resolution: " << resolution << std::endl;
}

/* Destructor */
FourierCurve::~FourierCurve()
{
    if (dev_coefficient_allocated)
    {
        cudaFree(dev_coefficients);
        std::cout << "dev_coefficients deallocated" << std::endl;
        dev_coefficient_allocated = false;
    }
    if (dev_trig_val_table_allocated)
    {
        cudaFree(dev_cos_table);
        cudaFree(dev_sin_table);
        std::cout << "dev_trig_val_table deallocated" << std::endl;
        dev_trig_val_table_allocated = false;
    }
    if (dev_curve_points_allocated)
    {
        cudaFree(dev_x);
        cudaFree(dev_y);
        cudaFree(dev_z);
        std::cout << "dev_curve_points deallocated" << std::endl;
        dev_curve_points_allocated = false;
    }
    if (dev_energy_matrix_allocated)
    {
        cudaFree(dev_energy_matrix);
        cudaFree(dev_differential_coefficients);
        std::cout << "dev_energy_matrix, dev_differential_coefficients deallocated" << std::endl;
        dev_energy_matrix_allocated = false;
    }

    std::cout << "Fourier curve destroyed!" << std::endl;
}

/* Load onto GPU */
void FourierCurve::cudafy()
{
    if (!dev_coefficient_allocated)
    {
        /* Allocate VRAM */
        cudaMalloc((void**) &dev_coefficients, sizeof(double) * 6 * (J + 1));

        dev_coefficient_allocated = true;
    }
    if (!dev_trig_val_table_allocated)
    {
        /* Allocate VRAM */
        cudaMalloc((void**) &dev_cos_table, sizeof(double) * (J + 1) * resolution);
        cudaMalloc((void**) &dev_sin_table, sizeof(double) * (J + 1) * resolution);

        /* Precomputed Trig Value Generation
           This table takes the form of

                                 | cos 0t | cos 1t | cos 2t | ... | cos Jt 
         t = 0                   |
         t = 2pi/res * 1         |
         t = 2pi/res * 2         |
         ...                     |
         t = 2pi/res * (res-1)   |

           But flattened.

           Access cos(k * 2pi/res * i) = dev_cos_table[k + (J + 1) i]
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

    /* Allocate scratch pad for x, y, z */
    if (!dev_curve_points_allocated)
    {
        cudaMalloc((void**) &dev_x, sizeof(double) * resolution);
        cudaMalloc((void**) &dev_y, sizeof(double) * resolution);
        cudaMalloc((void**) &dev_z, sizeof(double) * resolution);

        dev_curve_points_allocated = true;
    }

    if (!dev_energy_matrix_allocated)
    {
        cudaMalloc((void**) &dev_energy_matrix, sizeof(double) * resolution * resolution);
        /* dev_differential_coefficients is in the form of


            0   1   2   3   ...     6(J+1) - 1
        + | 
        - |

        Index by dev_differential_coefficients[6(J+1) * i + j]
        where i is 0 for plus and 1 for negative
        j is the coefficient index.
         */
        cudaMalloc((void**) &dev_differential_coefficients, sizeof(double) * 2 * 6 * (J + 1));

        dev_energy_matrix_allocated = true;
    }
                

    /* Copy all the coefficients, concatenated.
     At every iteration of cudafy function, this is the only part that is run again. */
    cudaMemcpy(dev_coefficients, &xa[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_coefficients + (J + 1), &xb[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_coefficients + 2 * (J + 1), &ya[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_coefficients + 3 * (J + 1), &yb[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_coefficients + 4 * (J + 1), &za[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_coefficients + 5 * (J + 1), &zb[0], sizeof(double) * (J + 1), cudaMemcpyHostToDevice);
}

void FourierCurve::cudaFlush()
{
    /* Curve points */
    cudaMemcpy(&x[0], dev_x, sizeof(double) * resolution, cudaMemcpyDeviceToHost);
    cudaMemcpy(&y[0], dev_y, sizeof(double) * resolution, cudaMemcpyDeviceToHost);
    cudaMemcpy(&z[0], dev_z, sizeof(double) * resolution, cudaMemcpyDeviceToHost);

    /* Coefficient values */
    cudaMemcpy(&xa[0], dev_coefficients, sizeof(double) * (J + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(&xb[0], dev_coefficients + (J + 1), sizeof(double) * (J + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ya[0], dev_coefficients + 2 * (J + 1), sizeof(double) * (J + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(&yb[0], dev_coefficients + 3 * (J + 1), sizeof(double) * (J + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(&za[0], dev_coefficients + 4 * (J + 1), sizeof(double) * (J + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(&zb[0], dev_coefficients + 5 * (J + 1), sizeof(double) * (J + 1), cudaMemcpyDeviceToHost);

    /* Differential */
    cudaMemcpy(&coeff_differential[0], dev_differential_coefficients, sizeof(double) * 6 * (J + 1), cudaMemcpyDeviceToHost);

    /* EnergyMatrix */
    cudaMemcpy(&energyMatrix[0], dev_energy_matrix, sizeof(double) * resolution * resolution, cudaMemcpyDeviceToHost);

    std::cout << "Flushed!" << std::endl;
}

__device__ void cross(double x1, double x2, double x3, double y1, double y2, double y3, double& z1, double& z2, double& z3)
{
    z1 = x2 * y3 - x3 * y2;
    z2 = x3 * y1 - x1 * y3;
    z3 = x1 * y2 - x2 * y1;
}

//Access cos(k * 2pi/res * i) = dev_cos_table[k + (J + 1) i]
__device__ double dev_trig_table_query(double* dev_table, unsigned int i, unsigned int k, unsigned int J)
{
    return dev_table[(J + 1) * i + k];
}

__device__ void fill_pos(double* dev_x, double* dev_y, double* dev_z, double* dev_coefficients, double* dev_cos_table, double* dev_sin_table, unsigned int resolution, unsigned int J)
{
    for (unsigned i = 0; i < resolution; i++)
    {
        dev_x[i] = dev_coefficients[0] / 2;
        dev_y[i] = dev_coefficients[(J + 1) * 2] / 2;
        dev_z[i] = dev_coefficients[(J + 1) * 4] / 2;
        for (unsigned k = 1; k <= J; k++)
        {
            dev_x[i] += dev_coefficients[k] * dev_trig_table_query(dev_cos_table, i, k, J);
            dev_x[i] += dev_coefficients[k + (J + 1)] * dev_trig_table_query(dev_sin_table, i, k, J);
            dev_y[i] += dev_coefficients[k + (J + 1) * 2] * dev_trig_table_query(dev_cos_table, i, k, J);
            dev_y[i] += dev_coefficients[k + (J + 1) * 3] * dev_trig_table_query(dev_sin_table, i, k, J);
            dev_z[i] += dev_coefficients[k + (J + 1) * 4] * dev_trig_table_query(dev_cos_table, i, k, J);
            dev_z[i] += dev_coefficients[k + (J + 1) * 5] * dev_trig_table_query(dev_sin_table, i, k, J);
        }
    }
}

__global__ void printCoefficientsPartiallyDEBUG(double* device_float_value)
{
    printf("%f\n", *device_float_value);
}

__global__ void crossDEBUG(double x1, double x2, double x3, double y1, double y2, double y3)
{
    double z1, z2, z3;
    cross(x1, x2, x3, y1, y2, y3, z1, z2, z3);
    printf("%f\n", z1);
    printf("%f\n", z2);
    printf("%f\n", z3);
}

__global__ void queryDEBUG(double* dev_table, int i, int k, unsigned int J)
{
    printf("%f\n", dev_trig_table_query(dev_table, i, k, J));
}

__global__ void fillDEBUG(double* dev_x, double* dev_y, double* dev_z, double* dev_coefficients, double* dev_cos_table, double* dev_sin_table, unsigned int resolution, unsigned int J)
{
    fill_pos(dev_x, dev_y, dev_z, dev_coefficients, dev_cos_table, dev_sin_table, resolution, J);
}
