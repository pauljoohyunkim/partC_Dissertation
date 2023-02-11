#ifndef __CURVEREPULSION_HPP__
#define __CURVEREPULSION_HPP__

#define ALPHA 2
#define BETA 4

#include "curve.hpp"

__device__ double kernelalphabeta(double px, double py, double pz, double qx, double qy, double qz, double Tx, double Ty, double Tz, double alpha = ALPHA, double beta = BETA);

__device__ double quadrature4PointSummand(double xix, double xiy, double xiz, double xipx, double xipy, double xipz, double xjx, double xjy, double xjz, double xjpx, double xjpy, double xjpz);

__global__ void tangentPointEnergyMatrixFill(double* dev_x, double* dev_y, double* dev_z, double* dev_energy_matrix, unsigned int resolution);

__global__ void sumEnergyMatrix(double* dev_energy_matrix, unsigned int resolution, double* dev_energy);

void fillDifferentialMatrix(FourierCurve& curve, double perturbation);

#endif
