#include "solver.hpp"
//#include "geometric-objects.hpp"

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
