#include "curve-repulsion-1.hpp"
#include "../solver.hpp"
#include "../geometric-objects.hpp"

#define DELTA_X 0.1
#define DELTA_T 0.05
#define LAMBDA 0.1

__global__ static void repulsiveCurveDifferential(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J);
__global__ static void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J);

int main()
{
    /* Generate curve */
    std::vector<double> x = { 1, 2, 3, 4, 5, 6 };
    std::vector<double> y = { 0, 2, 4, 6, 8, -1 };
    std::vector<double> z = { -1, -2, 3, -4, 5, 0.6 };
    cuRepulsiveCurve C(x, y, z);
    C.cudafy();


    /* Parallel Pool */
    dim3 grid(C.J, 3);

    /* Gradient Flow */
    repulsiveCurveDifferential<<<grid, 1>>>(C.dev_x, C.dev_y, C.dev_z, C.dev_energyMatrixFlattened, C.J);
    repulsiveCurveGradientFlow<<<grid, 1>>>(C.dev_x, C.dev_y, C.dev_z, C.dev_energyMatrixFlattened, C.J);
    C.flushFromDevice();

    return 0;
}

__global__ static void repulsiveCurveDifferential(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J)
{
    /* Point Index */
    int i = blockIdx.x;
    /* Coordinate Index j = 0, 1, 2 for x, y, z respectively */
    int j = blockIdx.y;
    //const double deltaX = 0.1;
    if (i < J)
    {
        unsigned int flattenedPos = J * j + i;
        if (j == 0)
        {
            cuDifferential(dev_x, dev_y, dev_z, i, DELTA_X, 0.0, 0.0, J, &dev_energyMatrixFlattened[flattenedPos]);
        }
        if (j == 1)
        {
            cuDifferential(dev_x, dev_y, dev_z, i, 0.0, DELTA_X, 0.0, J, &dev_energyMatrixFlattened[flattenedPos]);
        }
        if (j == 2)
        {
            cuDifferential(dev_x, dev_y, dev_z, i, 0.0, 0.0, DELTA_X, J, &dev_energyMatrixFlattened[flattenedPos]);
        }
    }
}

__global__ static void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J)
{
    /* Point Index */
    int i = blockIdx.x;
    /* Coordinate Index j = 0, 1, 2 for x, y, z respectively */
    int j = blockIdx.y;
    //const double deltaX = 0.1;
    if (i < J)
    {
        unsigned int flattenedPos = J * j + i;
        if (j == 0)
        {
            dev_x[i] += -dev_energyMatrixFlattened[flattenedPos] / DELTA_X * DELTA_T - LAMBDA * DELTA_T * dev_x[i];
        }
        if (j == 1)
        {
            dev_y[i] += -dev_energyMatrixFlattened[flattenedPos] / DELTA_X * DELTA_T - LAMBDA * DELTA_T * dev_y[i];
        }
        if (j == 2)
        {
            dev_z[i] += -dev_energyMatrixFlattened[flattenedPos] / DELTA_X * DELTA_T - LAMBDA * DELTA_T * dev_z[i];
        }
    }
}
