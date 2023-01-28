#include "curve-repulsion-1.hpp"
#include "../solver.hpp"
#include "../geometric-objects.hpp"

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
    repulsiveCurveGradientFlow<<<grid, 1>>>(C.dev_x, C.dev_y, C.dev_z, C.dev_energyMatrixFlattened, C.J);
    C.flushFromDevice();

    return 0;
}

__global__ static void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J)
{
    /* Point Index */
    int i = blockIdx.x;
    /* Coordinate Index j = 0, 1, 2 for x, y, z respectively */
    int j = blockIdx.y;
    //printf("repulsiveCurveGradientFlow(.): (i,j)=(%d,%d)\n", i, j);
    const double deltaX = 0.1;
    if (i < J)
    {
        unsigned int flattenedPos = J * j + i;
        if (j == 0)
        {
            dev_energyMatrixFlattened[flattenedPos] = cuDifferential(dev_x, dev_y, dev_z, i, deltaX, 0.0, 0.0, J);
        }
        if (j == 1)
        {
            dev_energyMatrixFlattened[flattenedPos] = cuDifferential(dev_x, dev_y, dev_z, i, 0.0, deltaX, 0.0, J);
        }
        if (j == 2)
        {
            dev_energyMatrixFlattened[flattenedPos] = cuDifferential(dev_x, dev_y, dev_z, i, 0.0, 0.0, deltaX, J);
        }
    }
}
