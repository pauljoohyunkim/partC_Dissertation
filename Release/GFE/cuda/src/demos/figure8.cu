#include <cmath>
#include <vector>
#include <iostream>
#include "figure8.cuh"
#include "../solver.cuh"
#include "../vector.cuh"
#include "../tpe.cuh"
#include "../export.cuh"

int main()
{
    /* Exporter */
    JsonExporter jsonX { "figure8_x.json" };
    JsonExporter jsonY { "figure8_y.json" };
    JsonExporter jsonZ { "figure8_z.json" };

    /* Figure 8 Generation */
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};
    for (auto i = 0; i < RESOLUTION; i++)
    {
        double theta = M_PI * 2 / RESOLUTION * i;
        x.push_back(SCALE * cos(theta));
        y.push_back(SCALE * sin(2 * theta));
        z.push_back(SCALE * TWIST * sin(theta));
    }

    /* Curve */
    CurveTensor figure8 { x, y, z };
    CurveTensor dCurve { RESOLUTION };

    /* Derivative Index */
    ScratchPad<int> s { RESOLUTION, 8 * (RESOLUTION - 3) };

    /* Evolve */
    for (auto t = 0; t < NUM_OF_STEPS; t++)
    {
        /* Energy Derivative Computaton */
        cuDEnergy<<<RESOLUTION, 1>>>(figure8.dev_blocks, dCurve.dev_blocks, s.scratchpads, RESOLUTION, ALPHA, BETA);
        cudaDeviceSynchronize();
        /* Evolution */
        tensorSubtract(figure8, dCurve, DELTA_T);
        //cudaDeviceSynchronize();
        /* Export to RAM, then to file */
        tensorBlockFlush(figure8, x, y, z);

        jsonX << x;
        jsonY << y;
        jsonZ << z;

        if (t % STATUS_UPDATE_FREQUENCY == 0)
        {
            std::cout << "Progress: " << t << "/" << NUM_OF_STEPS << std::endl;
        }
    }
}
