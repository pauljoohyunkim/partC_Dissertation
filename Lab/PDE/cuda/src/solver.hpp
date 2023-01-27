#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "geometric-objects.hpp"

#define ALPHA 2
#define BETA 4

/* Repulsive Curve Class */
class cuRepulsiveCurve: public cuCurve
{
    public:
        /* Constructor */
        cuRepulsiveCurve(unsigned int aJ);
        cuRepulsiveCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ);

        /* Deconstructor */
        ~cuRepulsiveCurve();

        double* dev_energyMatrix {};

        /* In order to use device functions, one needs to turn the array into CUDA array by
         * invoking cudafy in the beginning */
        void cudafy();
        /* Flush progress to x, y, z from GPU */
        void flushFromDevice();

    protected:
        bool dev_energyMatrix_allocated { false };
        std::vector<double> energyMatrixFlattened { };

};

__global__ void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrix, double J);
__global__ void repulsiveCurveEnergy(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrix, double J);

/* This function is used to fill energyMatrix; a matrix of summand for the tangent point energy.
 * Sum the values in the energyMatrix for the full tangent point energy of the curve.
 * */
__device__ void fillEnergyMatrix(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrix, unsigned int J);

__device__ void fillEnergyMatrixDifferential(double* dev_x, double* dev_y, double* dev_z, int index, double diffx, double diffy, double diffz, double* dev_energyMatrix, unsigned int J);

/* Pass coordinates of p, q, T */
__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double aAlpha, double aBeta);

/* Pass xi, x_{i+1}, x_j, x_{j+1}, Ti */
__device__ double kernelFunction(double xix, double xiy, double xiz, double xipx, double xipy, double xipz, double xjx, double xjy, double xjz, double xjpx, double xjpy, double xjpz, double Tix, double Tiy, double Tiz);

/* Pass a curve */
__device__ void cross(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3);

__device__ double l2norm3D(double x1, double x2, double x3);

#endif  // solver.hpp
