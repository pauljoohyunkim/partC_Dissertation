#include "curve-repulsion-1.hpp"
#include "../solver.hpp"
#include "../geometric-objects.hpp"
#include "../export.hpp"
#include <cmath>
#include <fstream>

#define DELTA_X 0.0005
#define DELTA_T 0.00001
#define LAMBDA 1
#define M 1000000

#define PLOT_FREQUENCY 10
#define AZIMUTHAL_SPEED 0.5
#define ELEVATION 3

#define PI 3.14159265358979

__global__ static void repulsiveCurveDifferential(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J);
__global__ static void repulsiveCurveGradientFlow(double* dev_x, double* dev_y, double* dev_z, double* dev_energyMatrixFlattened, double J);

int main()
{
    /* Generate curve */

    /* Example 1: */
    //std::vector<double> x = { 1, 2, 3, 4, 5, 6 };
    //std::vector<double> y = { 0, 2, 4, 6, 8, -1 };
    //std::vector<double> z = { -1, -2, 3, -4, 5, 0.6 };
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};

    /* Example 2: Helix and Semicircle */
    
    //const int resolution { 30 };
    //for (auto i = 0; i < resolution; i++)
    //{
    //    double theta = 4 * PI * (double) i / resolution;
    //    x.push_back(cos(theta));
    //    y.push_back(sin(theta));
    //    z.push_back(theta / (2 * PI));
    //}
    //for (auto i = 1; i < resolution; i++)
    //{
    //    double theta = PI * (double) i / resolution;
    //    x.push_back(1);
    //    y.push_back(2 * sin(theta));
    //    z.push_back(1 + cos(theta));
    //}

    /* Example 3: Just a circle */
    const int resolution { 32 };
    for (auto i = 0; i < resolution; i++)
    {
        double theta = 2 * PI * (double) i / resolution;
        x.push_back(cos(theta));
        y.push_back(sin(theta));
        z.push_back(0);
    }

    /* Example 4: Curve defined by three fourier series */
    //const int resolution { 80 };
    //for (auto i = 0; i < resolution; i++)
    //{
    //    double theta = 2 * PI * (double) i / resolution;
    //    x.push_back(cos(theta) + 3 * cos(2 * theta) - cos(3 * theta) - 0.7 * cos(4 * theta));
    //    y.push_back(sin(theta) + 0.2 * cos(2 * theta) + 2 * sin(3 * theta) + 0.2 * cos(4 * theta));
    //    z.push_back(cos(theta) + 2 * cos(2 * theta) + 4 * cos(3 * theta) - 2.2 * cos(4 * theta));
    //}

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


            /* Plotting */
            //auto curvePlot = matplot::plot3(C.x, C.y, C.z);
            //curvePlot->line_width(5);
            //matplot::view(AZIMUTHAL_SPEED * t, ELEVATION);
            //matplot::xrange({-5, 5});
            //matplot::yrange({-5, 5});
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
