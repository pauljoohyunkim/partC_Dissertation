#include <iostream>
#include "geometric-objects.hpp"
#include "solver.hpp"

/* Constructor for cuRepulsiveCurve */ 
cuRepulsiveCurve::cuRepulsiveCurve(unsigned int aJ): cuCurve(aJ){}


cuRepulsiveCurve::cuRepulsiveCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ): cuCurve(aX, aY, aZ){}

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

__device__ double cuRepulsiveCurve::energy()
{
    
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
    cross(px, py, pz, qx, qy, qz, numx, numy, numz);
    double numerator = pow(l2norm3D(numx, numy, numz), alpha);
    double denominator = pow(l2norm3D(pmqx, pmqy, pmqz), beta);

    return numerator / denominator;
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
