#include "../solver.hpp"
#include "../geometric-objects.hpp"
#include "../export.hpp"
#include <cmath>
#include <fstream>

#define DELTA_X 0.1
#define DELTA_T 0.005
#define LAMBDA 0.0001
#define M 10000

#define PLOT_FREQUENCY 10
#define AZIMUTHAL_SPEED 0.5
#define ELEVATION 3

__global__ static void repulsiveCurveDifferential(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J);
__global__ static void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J);

int main()
{
    /* Generate curve */

    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};

    const int resolution { 160 };
    for (auto i = 0; i < resolution; i++)
    {
        double theta = 2 * M_PI * (double) i / resolution;
        x.push_back(cos(theta) + 3 * cos(2 * theta) - cos(3 * theta) - 0.7 * cos(4 * theta));
        y.push_back(sin(theta) + 0.2 * cos(2 * theta) + 2 * sin(3 * theta) + 0.2 * cos(4 * theta));
        z.push_back(cos(theta) + 2 * cos(2 * theta) + 4 * cos(3 * theta) - 2.2 * cos(4 * theta));
    }

    cuRepulsiveCurve C(x, y, z);
    C.cudafy();


    /* Parallel Pool */
    dim3 grid(C.J, 3);

    /* Export to json */
    std::ofstream jsonX("x.json");
    std::ofstream jsonY("y.json");
    std::ofstream jsonZ("z.json");
    jsonX << "[";
    jsonY << "[";
    jsonZ << "[";

    /* Gradient Flow */
    for (auto t = 0; t < M; t++)
    {
        repulsiveCurveDifferential<<<grid, 1>>>(C.dev_x, C.dev_y, C.dev_z, C.dev_energyMatrixFlattened, C.J);
        repulsiveCurveGradientFlow<<<grid, 1>>>(C.dev_x, C.dev_y, C.dev_z, C.dev_energyMatrixFlattened, C.J);
        C.flushFromDevice();
        if (t % PLOT_FREQUENCY == 0)
        {
            std::cout << "Progress: " << t << "/" << M << "(" << (float) t / M * 100 << "%)" << std::endl;
            /* Export to json */
            if (t != 0)
            {
                jsonX << ",\n";
                jsonY << ",\n";
                jsonZ << ",\n";
            }
            vectorParse(jsonX, C.x, C.J);
            vectorParse(jsonY, C.y, C.J);
            vectorParse(jsonZ, C.z, C.J);
        }
    }
    jsonX << "]";
    jsonY << "]";
    jsonZ << "]";
    jsonX.close();
    jsonY.close();
    jsonZ.close();

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
            cuDifferentialSimple(dev_x, dev_y, dev_z, i, DELTA_X, 0.0, 0.0, J, &dev_energyMatrixFlattened[flattenedPos]);
        }
        if (j == 1)
        {
            cuDifferentialSimple(dev_x, dev_y, dev_z, i, 0.0, DELTA_X, 0.0, J, &dev_energyMatrixFlattened[flattenedPos]);
        }
        if (j == 2)
        {
            cuDifferentialSimple(dev_x, dev_y, dev_z, i, 0.0, 0.0, DELTA_X, J, &dev_energyMatrixFlattened[flattenedPos]);
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
