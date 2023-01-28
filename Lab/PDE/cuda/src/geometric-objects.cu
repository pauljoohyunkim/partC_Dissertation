#include "geometric-objects.hpp"
#include <iostream>
#include <exception>

/* Constructor for cuCurve */ 
cuCurve::cuCurve(unsigned int aJ)
{
    J = aJ;
    x.reserve(J);
    y.reserve(J);
    z.reserve(J);
    for (unsigned int i = 0; i < J; i++)
    {
        x.push_back(0.0);
        y.push_back(0.0);
        z.push_back(0.0);
    }
}

cuCurve::cuCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ)
{
    J = aX.size();
    x.reserve(J);
    y.reserve(J);
    z.reserve(J);
    for (unsigned int i = 0; i < J; i++)
    {
        x.push_back(aX[i]);
        y.push_back(aY[i]);
        z.push_back(aZ[i]);
    }
}

/* Deconstructor */
cuCurve::~cuCurve()
{
    if(dev_x_allocated)
    {
        cudaFree(dev_x);
    }
    if(dev_y_allocated)
    {
        cudaFree(dev_y);
    }
    if(dev_z_allocated)
    {
        cudaFree(dev_z);
    }

    std::cout << "cuCurve Deallocated" << std::endl;
}

/* Call this function before doing GPU stuff */
void cuCurve::cudafy()
{
    /* Allocate memory and copy data */
    cudaMalloc((void**)&dev_x, J * sizeof(double));
    dev_x_allocated = true;
    cudaMalloc((void**)&dev_y, J * sizeof(double));
    dev_y_allocated = true;
    cudaMalloc((void**)&dev_z, J * sizeof(double));
    dev_z_allocated = true;

    cudaMemcpy(dev_x, &x[0], J * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y, &y[0], J * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_z, &z[0], J * sizeof(double), cudaMemcpyHostToDevice);

    std::cout << "cuCurve Allocated" << std::endl;
}

void cuCurve::flushFromDevice()
{
    if (!dev_x_allocated)
    {
        throw std::runtime_error("cuCurve::flushFromDevice(): dev_x not allocated");
    }
    if (!dev_y_allocated)
    {
        throw std::runtime_error("cuCurve::flushFromDevice(): dev_y not allocated");
    }
    if (!dev_z_allocated)
    {
        throw std::runtime_error("cuCurve::flushFromDevice(): dev_z not allocated");
    }
    cudaMemcpy(&x[0], dev_x, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&y[0], dev_y, J * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&z[0], dev_z, J * sizeof(double), cudaMemcpyDeviceToHost);
}

__device__ double& cuCurve::getX(int i)
{
    return dev_x[((i % (int) J) + J) % J];
}

__device__ double& cuCurve::getY(int i)
{
    return dev_y[((i % (int) J) + J) % J];
}

__device__ double& cuCurve::getZ(int i)
{
    return dev_z[((i % (int) J) + J) % J];
}

/* coordnum = 0, 1, 2 for X, Y, Z respectively */
double cuCurve::getValFromDevice(int coordnum, int i)
{
    double val;
    i = ((i % (int) J) + J) % J;
    switch (coordnum)
    {
        /* X */
        case 0:
            cudaMemcpy(&val, &dev_x[i], sizeof(double), cudaMemcpyDeviceToHost);
            break;
        /* Y */
        case 1:
            cudaMemcpy(&val, &dev_y[i], sizeof(double), cudaMemcpyDeviceToHost);
            break;
        /* Z */
        case 2:
            cudaMemcpy(&val, &dev_z[i], sizeof(double), cudaMemcpyDeviceToHost);
            break;
        default:
            throw std::runtime_error("Coordnum out of bound.");
            break;
    }

    return val;
}
