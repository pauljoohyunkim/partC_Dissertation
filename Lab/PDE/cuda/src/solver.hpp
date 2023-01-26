#ifndef SOLVER_HPP
#define SOLVER_HPP


/* pass coordinates of p, q, T */
__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double aAlpha, double aBeta);
        
__device__ void cross(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3);
__device__ double l2norm3D(double x1, double x2, double x3);

#endif  // solver.hpp
