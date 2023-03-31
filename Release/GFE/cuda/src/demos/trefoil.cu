#include <cmath>
#include <vector>
#include <iostream>
#include "trefoil.cuh"
#include "../solver.cuh"
#include "../vector.cuh"
#include "../tpe.cuh"
#include "../export.cuh"

/* Filename Construction
   i = 0 for x
   i = 1 for y
   i = 2 for z
 */
static std::string filename(int i);

int main()
{
    /* Exporter */
    JsonExporter jsonX { filename(0) };
    JsonExporter jsonY { filename(1) };
    JsonExporter jsonZ { filename(2) };

    /* Figure 8 Generation */
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};
    for (auto i = 0; i < RESOLUTION; i++)
    {
        double theta = M_PI * 2 / RESOLUTION * i;
        x.push_back((sin(theta) + 2 * sin(2 * theta)));
        y.push_back((cos(theta) - 2 * cos(2 * theta)));
        z.push_back(-(1.2 - sin(theta)) * sin(3 * theta));
    }

    /* Curve */
    CurveTensor trefoil { x, y, z };
    CurveTensor dCurve { RESOLUTION };

    /* Derivative Index */
    ScratchPad<int> s { RESOLUTION, 8 * (RESOLUTION - 3) };

    /* Evolve */
    for (auto t = 0; t < NUM_OF_STEPS; t++)
    {
        /* Energy Derivative Computaton */
        cuDEnergy<<<RESOLUTION, 1>>>(trefoil.dev_blocks, dCurve.dev_blocks, s.scratchpads, RESOLUTION, ALPHA, BETA);
        cudaDeviceSynchronize();
        /* Evolution */
        tensorSubtract(trefoil, dCurve, DELTA_T);
        //cudaDeviceSynchronize();
        /* Export to RAM, then to file */
        tensorBlockFlush(trefoil, x, y, z);

        if (isnan(x[0]) || isnan(y[0]) || isnan(z[0]))
        {
            std::cout << "NaN Detected! Stopping now!" << std::endl;
            break;
        }

        if (t % STATUS_UPDATE_FREQUENCY == 0)
        {

            std::cout << "\33[2K\rProgress: " << t << "/" << NUM_OF_STEPS << " ";
            std::cout << "Tensor[0]: (" << x[0] << ", " << y[0] << ", " << z[0] << ")" << std::flush;
        }

        if (t % PLOT_FREQUENCY == 0)
        {
            jsonX << x;
            jsonY << y;
            jsonZ << z;
        }
    }
}

static std::string filename(int i)
{
    std::string fn { "trefoil" };
    fn += std::string("_ALPHA_") + std::to_string(ALPHA);
    fn += std::string("_BETA_") + std::to_string(BETA);
    fn += std::string("_RES_") + std::to_string(RESOLUTION);
    fn += std::string("_DELTAT_") + std::to_string(DELTA_T);
    fn += std::string("_NUMOFSTEP_") + std::to_string(NUM_OF_STEPS);
    fn += std::string("_PLOTFREQ_") + std::to_string(PLOT_FREQUENCY);
    fn += std::string("_SCALE_") + std::to_string(SCALE);
    fn += std::string("_TWIST_") + std::to_string(TWIST);
    if (i == 0)
    {
        fn += std::string("_x.json");
    }
    else if (i == 1)
    {
        fn += std::string("_y.json");
    }
    else if (i == 2)
    {
        fn += std::string("_z.json");
    }
    return fn;
}