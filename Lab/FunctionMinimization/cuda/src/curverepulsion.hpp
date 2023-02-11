#ifndef __CURVEREPULSION_HPP__
#define __CURVEREPULSION_HPP__

#define ALPHA 2
#define BETA 4

__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double alpha = ALPHA, double beta = BETA);

__device__ double quadrature4PointSummand(double xix, double xiy, double xiz, double xipx, double xipy, double xipz, double xjx, double xjy, double xjz, double xjpx, double xjpy, double xjpz);

__device__ double tangentPointEnergy(double* dev_x, double* dev_y, double* dev_z, unsigned int resolution);

__global__ void energyDEBUG(double* dev_x, double* dev_y, double* dev_z, unsigned int resolution);

#endif
