#include "curverepulsion.hpp"
#include "curve.hpp"
#include <iostream>

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

//__global__ void energyDEBUG(double* dev_x, double* dev_y, double* dev_z, unsigned int resolution)
//{
//    double energy = tangentPointEnergy(dev_x, dev_y, dev_z, resolution);
//    printf("Energy: %f\n", energy);
//}
