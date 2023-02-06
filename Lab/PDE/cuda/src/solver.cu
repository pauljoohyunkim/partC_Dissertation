#include <iostream>
#include "geometric-objects.hpp"
#include "solver.hpp"

/* Constructor for cuRepulsiveCurve */ 
cuRepulsiveCurve::cuRepulsiveCurve(unsigned int aJ): cuCurve(aJ)
{
    energyMatrixFlattened.resize(3 * aJ);
}


cuRepulsiveCurve::cuRepulsiveCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ): cuCurve(aX, aY, aZ)
{
    unsigned int aJ = aX.size();
    energyMatrixFlattened.resize(3 * aJ);
}

/* Deconstructor */
cuRepulsiveCurve::~cuRepulsiveCurve()
{
    if (dev_x_allocated)
    {
        cudaFree(dev_x);
        dev_x_allocated = false;
    }
    if (dev_y_allocated)
    {
        cudaFree(dev_y);
        dev_y_allocated = false;
    }
    if (dev_z_allocated)
    {
        cudaFree(dev_z);
        dev_z_allocated = false;
    }
    if (dev_energyMatrixFlattened_allocated)
    {
        cudaFree(dev_energyMatrixFlattened);
        dev_energyMatrixFlattened_allocated = false;
    }

    std::cout << "cuRepulsiveCurve Deallocated" << std::endl;
}

/* Call this function before doing GPU stuff */
void cuRepulsiveCurve::cudafy()
{
    /* Allocate memory */
    cudaMalloc((void**)&dev_x, J * sizeof(double));
    dev_x_allocated = true;
    cudaMalloc((void**)&dev_y, J * sizeof(double));
    dev_y_allocated = true;
    cudaMalloc((void**)&dev_z, J * sizeof(double));
    dev_z_allocated = true;
    cudaMalloc((void**)&dev_energyMatrixFlattened, 3 * J * sizeof(double));
    dev_energyMatrixFlattened_allocated = true;

    /* Copy curve points */
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
    if (!dev_energyMatrixFlattened_allocated)
    {
        throw std::runtime_error("cuRepulsiveCurve::flushFromDevice: dev_energyMatrixFlattened not allocated");
    }
    cudaMemcpy(&x[0], dev_x, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&y[0], dev_y, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&z[0], dev_z, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&energyMatrixFlattened[0], dev_energyMatrixFlattened, 3 * J * sizeof(double), cudaMemcpyDeviceToHost);
    //cudaMemcpy(&energyMatrixFlattened[0], &dev_energyTensor_x[0], J * J * sizeof(double), cudaMemcpyDeviceToHost);
    //cudaMemcpy(&energyMatrixFlattened[0], dev_energyMatrix, J * J * sizeof(double), cudaMemcpyDeviceToHost);
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
    //printf("fillEnergyMatrix() called.\n");
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
            //printf("i: %d, j: %d, energyLocal = %f\n", i, j, dev_energyMatrix[flattenPos]);
        }
        else
        {
            dev_energyMatrix[flattenPos] = 0;
        }
        }
    }
    //}
}

__device__ double cuDifferential(double* dev_x, double* dev_y, double* dev_z, int index, double diffx, double diffy, double diffz, unsigned int J)
{
    /* 2-Pass:
     On the first pass, it perturbs the curve a bit temporarily and computes the energy.
     On the second pass, it perturbs the curve in the opposite direction, then computes the energy, subtracting off kernel points
     then it divides by 2.*/

    double differential { 0 };

    index = ((index % (int) J) + J) % J;

    /* Energy of perturbed curve */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
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
                /* Perturbation */
                if (i == index)
                {
                    xix += diffx;
                    xiy += diffy;
                    xiz += diffz;
                }
                if (j == index)
                {
                    xjx += diffx;
                    xjy += diffy;
                    xjz += diffz;
                }

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];
                
                /* Perturbation */
                if (ip1 == index)
                {
                    xipx += diffx;
                    xipy += diffy;
                    xipz += diffz;
                }
                if (jp1 == index)
                {
                    xjpx += diffx;
                    xjpy += diffy;
                    xjpz += diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                /* TI = pI / lI */
                double TIx = xIx / lI;
                double TIy = xIy / lI;
                double TIz = xIz / lI;

                differential += kernelFunction(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz, TIx, TIy, TIz) * lI * lJ;
            }
        }
    }

    /* Energy of curve perturbed in the opposite direction subtracted off */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
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
                /* Perturbation */
                if (i == index)
                {
                    xix -= diffx;
                    xiy -= diffy;
                    xiz -= diffz;
                }
                if (j == index)
                {
                    xjx -= diffx;
                    xjy -= diffy;
                    xjz -= diffz;
                }

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];
                /* Perturbation */
                if (ip1 == index)
                {
                    xipx -= diffx;
                    xipy -= diffy;
                    xipz -= diffz;
                }
                if (jp1 == index)
                {
                    xjpx -= diffx;
                    xjpy -= diffy;
                    xjpz -= diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                /* TI = pI / lI */
                double TIx = xIx / lI;
                double TIy = xIy / lI;
                double TIz = xIz / lI;

                differential -= kernelFunction(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz, TIx, TIy, TIz) * lI * lJ;
            }
        }
    }
    
    return differential / 2;
}


__device__ void cuDifferential(double* dev_x, double* dev_y, double* dev_z, int index, double diffx, double diffy, double diffz, unsigned int J, double* pVar)
{
    /* 2-Pass:
     On the first pass, it perturbs the curve a bit temporarily and computes the energy.
     On the second pass, it perturbs the curve in the opposite direction, then computes the energy, subtracting off kernel points
     then it divides by 2.*/

    *pVar = 0;

    index = ((index % (int) J) + J) % J;

    /* Energy of perturbed curve */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
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
                /* Perturbation */
                if (i == index)
                {
                    xix += diffx;
                    xiy += diffy;
                    xiz += diffz;
                }
                if (j == index)
                {
                    xjx += diffx;
                    xjy += diffy;
                    xjz += diffz;
                }

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];
                
                /* Perturbation */
                if (ip1 == index)
                {
                    xipx += diffx;
                    xipy += diffy;
                    xipz += diffz;
                }
                if (jp1 == index)
                {
                    xjpx += diffx;
                    xjpy += diffy;
                    xjpz += diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                /* TI = pI / lI */
                double TIx = xIx / lI;
                double TIy = xIy / lI;
                double TIz = xIz / lI;

                *pVar += kernelFunction(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz, TIx, TIy, TIz) * lI * lJ;
            }
        }
    }

    /* Energy of curve perturbed in the opposite direction subtracted off */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
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
                /* Perturbation */
                if (i == index)
                {
                    xix -= diffx;
                    xiy -= diffy;
                    xiz -= diffz;
                }
                if (j == index)
                {
                    xjx -= diffx;
                    xjy -= diffy;
                    xjz -= diffz;
                }

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];
                /* Perturbation */
                if (ip1 == index)
                {
                    xipx -= diffx;
                    xipy -= diffy;
                    xipz -= diffz;
                }
                if (jp1 == index)
                {
                    xjpx -= diffx;
                    xjpy -= diffy;
                    xjpz -= diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                /* TI = pI / lI */
                double TIx = xIx / lI;
                double TIy = xIy / lI;
                double TIz = xIz / lI;

                *pVar -= kernelFunction(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz, TIx, TIy, TIz) * lI * lJ;
            }
        }
    }
    *pVar = (*pVar) / 2;
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


__device__ double cuDifferentialSimple(double* dev_x, double* dev_y, double* dev_z, int index, double diffx, double diffy, double diffz, unsigned int J, double* pVar)
{
    double de { 0 };
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
            if (i != j)
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

                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];

                /* Perturbation */
                if (i == index)
                {
                    xix += diffx;
                    xiy += diffy;
                    xiz += diffz;
                }
                if (j == index)
                {
                    xjx += diffx;
                    xjy += diffy;
                    xjz += diffz;
                }
                if (ip1 == index)
                {
                    xipx += diffx;
                    xipy += diffy;
                    xipz += diffz;
                }
                if (jp1 == index)
                {
                    xjpx += diffx;
                    xjpy += diffy;
                    xjpz += diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                de += lI * lJ / l2norm3D(xix - xjx , xiy - xjy, xiz - xjz);
            }
        }
    }


    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
            if (i != j)
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
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];

                /* Perturbation */
                if (i == index)
                {
                    xix -= diffx;
                    xiy -= diffy;
                    xiz -= diffz;
                }
                if (j == index)
                {
                    xjx -= diffx;
                    xjy -= diffy;
                    xjz -= diffz;
                }
                if (ip1 == index)
                {
                    xipx -= diffx;
                    xipy -= diffy;
                    xipz -= diffz;
                }
                if (jp1 == index)
                {
                    xjpx -= diffx;
                    xjpy -= diffy;
                    xjpz -= diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                de -= lI * lJ / l2norm3D(xix - xjx , xiy - xjy, xiz - xjz);
            }
        }
    }

    *pVar = de / 2;
}



__device__ void cuDifferentialSimpleCentral(double* dev_x, double* dev_y, double* dev_z, int index, double diffx, double diffy, double diffz, unsigned int J, double* pVar)
{
    /* 2-Pass:
     On the first pass, it perturbs the curve a bit temporarily and computes the energy.
     On the second pass, it perturbs the curve in the opposite direction, then computes the energy, subtracting off kernel points
     then it divides by 2.*/

    *pVar = 0;

    index = ((index % (int) J) + J) % J;

    /* Energy of perturbed curve */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
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
                /* Perturbation */
                if (i == index)
                {
                    xix += diffx;
                    xiy += diffy;
                    xiz += diffz;
                }
                if (j == index)
                {
                    xjx += diffx;
                    xjy += diffy;
                    xjz += diffz;
                }

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];
                
                /* Perturbation */
                if (ip1 == index)
                {
                    xipx += diffx;
                    xipy += diffy;
                    xipz += diffz;
                }
                if (jp1 == index)
                {
                    xjpx += diffx;
                    xjpy += diffy;
                    xjpz += diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                *pVar += kernelFunctionSimpleCentral(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz) * lI * lJ;
            }
        }
    }

    /* Energy of curve perturbed in the opposite direction subtracted off */
    for (int i = 0; i < J; i++)
    {
        for (int j = 0; j < J; j++)
        {
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
                /* Perturbation */
                if (i == index)
                {
                    xix -= diffx;
                    xiy -= diffy;
                    xiz -= diffz;
                }
                if (j == index)
                {
                    xjx -= diffx;
                    xjy -= diffy;
                    xjz -= diffz;
                }

                /* x_{i+1}, x_{j+1} */
                double xipx = dev_x[ip1];
                double xipy = dev_y[ip1];
                double xipz = dev_z[ip1];
                double xjpx = dev_x[jp1];
                double xjpy = dev_y[jp1];
                double xjpz = dev_z[jp1];
                /* Perturbation */
                if (ip1 == index)
                {
                    xipx -= diffx;
                    xipy -= diffy;
                    xipz -= diffz;
                }
                if (jp1 == index)
                {
                    xjpx -= diffx;
                    xjpy -= diffy;
                    xjpz -= diffz;
                }

                /* xI, xJ */
                double xIx = xipx - xix;
                double xIy = xipy - xiy;
                double xIz = xipz - xiz;
                double xJx = xjpx - xjx;
                double xJy = xjpy - xjy;
                double xJz = xjpz - xjz;

                /* lI, lJ */
                double lI = l2norm3D(xIx, xIy, xIz);
                double lJ = l2norm3D(xJx, xJy, xJz);

                *pVar -= kernelFunctionSimpleCentral(xix, xiy, xiz, xipx, xipy, xipz,
                        xjx, xjy, xjz, xjpx, xjpy, xjpz) * lI * lJ;
            }
        }
    }
    *pVar = (*pVar) / 2;
}


__device__ double kernelFunctionSimpleCentral(double xix, double xiy, double xiz, double xipx, double xipy, double xipz, double xjx, double xjy, double xjz, double xjpx, double xjpy, double xjpz, double alpha)
{
    double kernelVal { 0 };
    kernelVal += 1 / pow(l2norm3D(xix - xjx, xiy - xjy, xiz - xjz), alpha);
    kernelVal += 1 / pow(l2norm3D(xix - xjpx, xiy - xjpy, xiz - xjpz), alpha);
    kernelVal += 1 / pow(l2norm3D(xipx - xjx, xipy - xjy, xipz - xjz), alpha);
    kernelVal += 1 / pow(l2norm3D(xipx - xjpx, xipy - xjpy, xipz - xjpz), alpha);

    return kernelVal / 4;
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
