#include "curverepulsion.hpp"
#include "curve.hpp"
#include <iostream>

__global__ static void centralDifference(double& energyplus, double& energyminus, double perturbation);

__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double alpha, double beta)
{
    double pmqx = px - qx;
    double pmqy = py - qy;
    double pmqz = pz - qz;

    double nx, ny, nz;
    cross(Tx, Ty, Tz, pmqx, pmqy, pmqz, nx, ny, nz);

    double numerator = norm3d(nx, ny, nz);
    double denominator = norm3d(pmqx, pmqy, pmqz);

    return pow(numerator, alpha) / pow(denominator, beta);
}

__device__ double quadrature4PointSummand(double xix, double xiy, double xiz, double xipx, double xipy, double xipz, double xjx, double xjy, double xjz, double xjpx, double xjpy, double xjpz)
{
    double xIx = xipx - xix;
    double xIy = xipy - xiy;
    double xIz = xipz - xiz;
    double xINorm = norm3d(xIx, xIy, xIz);

    double xJx = xjpx - xjx;
    double xJy = xjpy - xjy;
    double xJz = xjpz - xjz;
    double xJNorm = norm3d(xJx, xJy, xJz);

    double TIx = xIx / xINorm;
    double TIy = xIy / xINorm;
    double TIz = xIz / xINorm;

    double res { 0 };

    res += kernelalphabeta(xix, xiy, xiz, xjx, xjy, xjz, TIx, TIy, TIz);
    res += kernelalphabeta(xix, xiy, xiz, xjpx, xjpy, xjpz, TIx, TIy, TIz);
    res += kernelalphabeta(xipx, xipy, xipz, xjx, xjy, xjz, TIx, TIy, TIz);
    res += kernelalphabeta(xipx, xipy, xipz, xjpx, xjpy, xjpz, TIx, TIy, TIz);

    return res / 4 * xINorm * xJNorm;
}

/* Pass in resolution x resolution */
__global__ void tangentPointEnergyMatrixFill(double* dev_x, double* dev_y, double* dev_z, double* dev_energy_matrix, unsigned int resolution)
{
    int i = blockIdx.x;
    int j = blockIdx.y;
    if (i < (int) resolution)
    {
        if (j < (int) resolution)
        {
            if (abs(i - j) > 1 && abs(i - j + (int) resolution) > 1 && abs(i - j - (int) resolution) > 1)
            {
                /* i+1, j+1 */
                int ip = (i + 1) % resolution;
                int jp = (j + 1) % resolution;

                /* x_i, x_j */
                double xix = dev_x[i];
                double xiy = dev_y[i];
                double xiz = dev_z[i];
                double xjx = dev_x[j];
                double xjy = dev_y[j];
                double xjz = dev_z[j];

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip];
                double xipy = dev_y[ip];
                double xipz = dev_z[ip];
                double xjpx = dev_x[jp];
                double xjpy = dev_y[jp];
                double xjpz = dev_z[jp];

                dev_energy_matrix[resolution * i + j] = quadrature4PointSummand(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz);
            }
        }
    }
}

/* Pass in 1 for parallel blocks */
__global__ void sumEnergyMatrix(double* dev_energy_matrix, unsigned int resolution, double* dev_energy)
{
    double energy { 0 };
    for (unsigned int i = 0; i < resolution * resolution; i++)
    {
        energy += dev_energy_matrix[i];
    }

    *dev_energy = energy;
}

__global__ static void centralDifference(double& energyplus, double& energyminus, double perturbation)
{
    energyplus = energyplus - energyminus;
    energyplus = energyplus / (2 * perturbation);
}


/* After this operation, first out of the two rows of curve.dev_differential_coefficients will be gradient of tangent point energy */
void fillDifferentialMatrix(FourierCurve& curve, double perturbation)
{
    /* xa perturbation */
    for (unsigned int coeffIndex = 0; coeffIndex < 6 * (curve.J + 1); coeffIndex++)
    {
        double temp;
        temp = curve.xa[coeffIndex];
        /* Perturbation+ */
        curve.xa[coeffIndex] += perturbation;
        curve.cudafy();
        fill_pos_from_host<<<1,1>>>(curve.dev_x, curve.dev_y, curve.dev_z, curve.dev_coefficients, curve.dev_cos_table, curve.dev_sin_table, curve.resolution, curve.J);
        dim3 grid(curve.resolution, curve.resolution);
        tangentPointEnergyMatrixFill<<<grid, 1>>>(curve.dev_x, curve.dev_y, curve.dev_z, curve.dev_energy_matrix, curve.resolution);
        sumEnergyMatrix<<<1,1>>>(curve.dev_energy_matrix, curve.resolution, curve.dev_differential_coefficients + coeffIndex);
        //printCoefficientsPartiallyDEBUG<<<1,1>>>(&curve.dev_differential_coefficients[coeffIndex]);
        cudaDeviceSynchronize();

        curve.xa[coeffIndex] = temp;

        /* Perturbation- */
        curve.xa[coeffIndex] -= perturbation;
        curve.cudafy();
        fill_pos_from_host<<<1,1>>>(curve.dev_x, curve.dev_y, curve.dev_z, curve.dev_coefficients, curve.dev_cos_table, curve.dev_sin_table, curve.resolution, curve.J);
        tangentPointEnergyMatrixFill<<<grid, 1>>>(curve.dev_x, curve.dev_y, curve.dev_z, curve.dev_energy_matrix, curve.resolution);
        sumEnergyMatrix<<<1,1>>>(curve.dev_energy_matrix, curve.resolution, curve.dev_differential_coefficients + 6 * (curve.J + 1) + coeffIndex);
        //printCoefficientsPartiallyDEBUG<<<1,1>>>(&curve.dev_differential_coefficients[6 * (curve.J + 1) + coeffIndex]);
        cudaDeviceSynchronize();

        curve.xa[coeffIndex] = temp;
        curve.cudafy();

        centralDifference<<<1,1>>>(curve.dev_differential_coefficients[coeffIndex], curve.dev_differential_coefficients[6 * (curve.J + 1) + coeffIndex], perturbation);

        //printCoefficientsPartiallyDEBUG<<<1,1>>>(&curve.dev_differential_coefficients[coeffIndex]);
        cudaDeviceSynchronize();

        printf("Done\n");
    }
}

/* Pass in 6 (J + 1) */
__global__ void gradientDescent(double* dev_coefficients, double* dev_differential_coefficients, double stepsize, unsigned int J)
{
    int i = blockIdx.x;
    if (i < 6 * ((int) J + 1))
    {
        dev_coefficients[i] -= stepsize * dev_differential_coefficients[i];
    }
    
}
