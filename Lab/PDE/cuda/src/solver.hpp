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

        __device__ double energy();

        /* In order to use device functions, one needs to turn the array into CUDA array by
         * invoking cudafy in the beginning */
        void cudafy();

    protected:
        double* dev_energyMatrix {};
        bool dev_energyMatrix_allocated { false };

};

/* Pass coordinates of p, q, T */
__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double aAlpha, double aBeta);

/* Pass a curve */
__device__ void cross(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3);

__device__ double l2norm3D(double x1, double x2, double x3);

#endif  // solver.hpp
