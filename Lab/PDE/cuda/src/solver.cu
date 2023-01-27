#include <iostream>
#include "geometric-objects.hpp"
#include "solver.hpp"

/* Constructor for cuRepulsiveCurve */ 
cuRepulsiveCurve::cuRepulsiveCurve(unsigned int aJ): cuCurve(aJ)
{
    energyMatrixFlattened.resize(aJ * aJ);
}


cuRepulsiveCurve::cuRepulsiveCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ): cuCurve(aX, aY, aZ)
{
    double aJ = aX.size();
    energyMatrixFlattened.resize(aJ * aJ);
}

/* Deconstructor */
cuRepulsiveCurve::~cuRepulsiveCurve()
{
    if (dev_x_allocated)
    {
        cudaFree(dev_x);
    }
    if (dev_y_allocated)
    {
        cudaFree(dev_y);
    }
    if (dev_z_allocated)
    {
        cudaFree(dev_z);
    }
    if (dev_energyMatrix_allocated)
    {
        cudaFree(dev_energyMatrix);
    }

    std::cout << "cuRepulsiveCurve Deallocated" << std::endl;
}

/* Call this function before doing GPU stuff */
void cuRepulsiveCurve::cudafy()
{
    /* Allocate memory and copy data */
    cudaMalloc((void**)&dev_x, J * sizeof(double));
    dev_x_allocated = true;
    cudaMalloc((void**)&dev_y, J * sizeof(double));
    dev_y_allocated = true;
    cudaMalloc((void**)&dev_z, J * sizeof(double));
    dev_z_allocated = true;
    cudaMalloc((void**)&dev_energyMatrix, (J * J) * sizeof(double*));
    dev_energyMatrix_allocated = true;

    cudaMemcpy(dev_x, &x[0], J * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y, &y[0], J * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_z, &z[0], J * sizeof(double), cudaMemcpyHostToDevice);

    std::cout << "cuRepulsiveCurve Allocated" << std::endl;
}

void cuRepulsiveCurve::flushFromDevice()
{
    if (!dev_x_allocated)
    {
        throw std::runtime_error("cuRepulsiveCurve::flushFromDevice: dev_x not allocated");
    }
    if (!dev_y_allocated)
    {
        throw std::runtime_error("cuRepulsiveCurve::flushFromDevice: dev_y not allocated");
    }
    if (!dev_z_allocated)
    {
        throw std::runtime_error("cuRepulsiveCurve::flushFromDevice: dev_z not allocated");
    }
    if (!dev_energyMatrix_allocated)
    {
        throw std::runtime_error("cuRepulsiveCurve::flushFromDevice: dev_energyMatrix not allocated");
    }
    cudaMemcpy(&x[0], dev_x, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&y[0], dev_y, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&z[0], dev_z, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&energyMatrixFlattened[0], dev_energyMatrix, J * J * sizeof(double), cudaMemcpyDeviceToHost);

}

__global__ void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrix, double J)
{
    fillEnergyMatrixDifferential(dev_x, dev_y, dev_z, 0, 0.1, 0.0, 0.0, dev_energyMatrix, J);
}

// Test Function
__global__ void repulsiveCurveEnergy(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrix, double J)
{
    fillEnergyMatrix(dev_x, dev_y, dev_z, dev_energyMatrix, J);
}

__device__ double sumArray(double* arr, unsigned int length)
{
    double sum = 0;
    for (unsigned int i = 0; i < length; i++)
    {
        sum += arr[i];
    }

    return sum;
}

__device__ void fillEnergyMatrix(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrix, unsigned int J)
{
    //int i = blockIdx.x;
    //int ip1 = (i + 1) % J;
    //int j = blockIdx.y;
    //int jp1 = (j + 1) % J;
    //if (i < (int) J && j < (int) J)
    //{
    printf("fillEnergyMatrix() called.\n");
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {

        int flattenPos = i + J * j;

        if (abs(i - j) > 1 && abs(i - j + (int) J) > 1 && abs(i - j - (int) J) > 1)
        {
            int ip1 = (i + 1) % J;
            int jp1 = (j + 1) % J;
            /* x_i, x_j */
            double xix = dev_x[i];
            double xiy = dev_y[i];
            double xiz = dev_z[i];
            double xjx = dev_x[j];
            double xjy = dev_y[j];
            double xjz = dev_z[j];

            /* xI, xJ */
            double xIx = dev_x[ip1] - xix;
            double xIy = dev_y[ip1] - xiy;
            double xIz = dev_z[ip1] - xiz;
            double xJx = dev_x[jp1] - xjx;
            double xJy = dev_y[jp1] - xjy;
            double xJz = dev_z[jp1] - xjz;

            /* lI, lJ */
            double lI = l2norm3D(xIx, xIy, xIz);
            double lJ = l2norm3D(xJx, xJy, xJz);

            /* TI = pI / lI */
            double TIx = xIx / lI;
            double TIy = xIy / lI;
            double TIz = xIz / lI;

            dev_energyMatrix[flattenPos] = kernelFunction(xix, xiy, xiz, dev_x[ip1], dev_y[ip1], dev_z[ip1],
                    xjx, xjy, xjz, dev_x[jp1], dev_y[jp1], dev_z[jp1], TIx, TIy, TIz) * lI * lJ;
            printf("i: %d, j: %d, energyLocal = %f\n", i, j, dev_energyMatrix[flattenPos]);
        }
        else
        {
            dev_energyMatrix[flattenPos] = 0;
        }
        }
    }
    //}
}

__device__ void fillEnergyMatrixDifferential(double* dev_x, double* dev_y, double* dev_z, int index, double diffx, double diffy, double diffz, double* dev_energyMatrix, unsigned int J)
{
    /* 2-Pass:
     On the first pass, it perturbs the curve a bit temporarily and computes the energy.
     On the second pass, it restores the original curve, then computes the energy, subtracting off kernel points*/

    printf("fillEnergyMatrixDifferential() called.\n");

    index = ((index % (int) J) + J) % J;

    double save_x = dev_x[index];
    double save_y = dev_y[index];
    double save_z = dev_z[index];


    /* Perturbation */
    dev_x[index] += diffx;
    dev_y[index] += diffy;
    dev_z[index] += diffz;

    /* Energy of perturbed curve */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {

            int flattenPos = i + J * j;

            if (abs(i - j) > 1 && abs(i - j + (int) J) > 1 && abs(i - j - (int) J) > 1)
            {
                int ip1 = (i + 1) % J;
                int jp1 = (j + 1) % J;
                /* x_i, x_j */
                double xix = dev_x[i];
                double xiy = dev_y[i];
                double xiz = dev_z[i];
                double xjx = dev_x[j];
                double xjy = dev_y[j];
                double xjz = dev_z[j];

                /* xI, xJ */
                double xIx = dev_x[ip1] - xix;
                double xIy = dev_y[ip1] - xiy;
                double xIz = dev_z[ip1] - xiz;
                double xJx = dev_x[jp1] - xjx;
                double xJy = dev_y[jp1] - xjy;
                double xJz = dev_z[jp1] - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                /* TI = pI / lI */
                double TIx = xIx / lI;
                double TIy = xIy / lI;
                double TIz = xIz / lI;

                dev_energyMatrix[flattenPos] = kernelFunction(xix, xiy, xiz, dev_x[ip1], dev_y[ip1], dev_z[ip1],
                        xjx, xjy, xjz, dev_x[jp1], dev_y[jp1], dev_z[jp1], TIx, TIy, TIz) * lI * lJ;
                printf("i: %d, j: %d, energyLocal = %f\n", i, j, dev_energyMatrix[flattenPos]);
            }
            else
            {
                dev_energyMatrix[flattenPos] = 0;
            }
        }
    }

    /* Restore Curve */
    dev_x[index] = save_x;
    dev_y[index] = save_y;
    dev_z[index] = save_z;
    /* Energy of original curve subtracted off */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {

            int flattenPos = i + J * j;

            if (abs(i - j) > 1 && abs(i - j + (int) J) > 1 && abs(i - j - (int) J) > 1)
            {
                int ip1 = (i + 1) % J;
                int jp1 = (j + 1) % J;
                /* x_i, x_j */
                double xix = dev_x[i];
                double xiy = dev_y[i];
                double xiz = dev_z[i];
                double xjx = dev_x[j];
                double xjy = dev_y[j];
                double xjz = dev_z[j];

                /* xI, xJ */
                double xIx = dev_x[ip1] - xix;
                double xIy = dev_y[ip1] - xiy;
                double xIz = dev_z[ip1] - xiz;
                double xJx = dev_x[jp1] - xjx;
                double xJy = dev_y[jp1] - xjy;
                double xJz = dev_z[jp1] - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                /* TI = pI / lI */
                double TIx = xIx / lI;
                double TIy = xIy / lI;
                double TIz = xIz / lI;

                dev_energyMatrix[flattenPos] -= kernelFunction(xix, xiy, xiz, dev_x[ip1], dev_y[ip1], dev_z[ip1],
                        xjx, xjy, xjz, dev_x[jp1], dev_y[jp1], dev_z[jp1], TIx, TIy, TIz) * lI * lJ;
                printf("i: %d, j: %d, energyLocalDifferential = %f\n", i, j, dev_energyMatrix[flattenPos]);
            }
            else
            {
                dev_energyMatrix[flattenPos] = 0;
            }
        }
    }
    
}

__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double alpha, double beta)
{
    double pmqx = px - qx;
    double pmqy = py - qy;
    double pmqz = pz - qz;
    double numx;
    double numy;
    double numz;

    /* T x (p-q) */
    cross(Tx, Ty, Tz, pmqx, pmqy, pmqz, numx, numy, numz);
    double numerator = pow(l2norm3D(numx, numy, numz), alpha);
    double denominator = pow(l2norm3D(pmqx, pmqy, pmqz), beta);

    return numerator / denominator;
}


__device__ double kernelFunction(double xix, double xiy, double xiz, double xipx, double xipy, double xipz, double xjx, double xjy, double xjz, double xjpx, double xjpy, double xjpz, double Tix, double Tiy, double Tiz)
{
    double kij { 0 };

    kij += kernelalphabeta(xix, xiy, xiz, xjx, xjy, xjz, Tix, Tiy, Tiz, ALPHA, BETA);
    kij += kernelalphabeta(xix, xiy, xiz, xjpx, xjpy, xjpz, Tix, Tiy, Tiz, ALPHA, BETA);
    kij += kernelalphabeta(xipx, xipy, xipz, xjx, xjy, xjz, Tix, Tiy, Tiz, ALPHA, BETA);
    kij += kernelalphabeta(xipx, xipy, xipz, xjpx, xjpy, xjpz, Tix, Tiy, Tiz, ALPHA, BETA);

    return kij / 4;
}

__device__ void cross(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3)
{
    x3 = y1 * z2 - y2 * z1;
    y3 = z1 * x2 - x1 * z2;
    z3 = x1 * y2 - x2 * y1;
}

__device__ double l2norm3D(double x1, double x2, double x3)
{
    double norm { 0 };
    norm += x1 * x1;
    norm += x2 * x2;
    norm += x3 * x3;

    return sqrt(norm);
}
