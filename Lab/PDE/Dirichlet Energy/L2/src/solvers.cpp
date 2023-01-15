#include "solvers.hpp"
#include "params.hpp"
#include <functional>
#include <iostream>
#include <vector>
#include <matplot/matplot.h>

void heatEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, unsigned int T, unsigned int M)
{
    double deltaT = (double) T / M;
    double deltaX = (b - a) / J;
    double cfl = deltaT / (deltaX * deltaX);
    double *u = new double [J + 1];
    std::vector<double> x {};
    double *ump1 = new double [J + 1];
    std::vector<double> plotVec(J + 1, 0);
    
    if (cfl > 0.5)
    {
        std::cout << "Warning: the CFL number is " << cfl << ", so stability is not guaranteed";
    }

    /* Initial Datum */
    for (auto j = 0; j <= J; j++)
    {
        x.push_back(a + (double) j * deltaX);
        u[j] = u_initial(x[j]);
    }


    matplot::hold(matplot::on);
    
    /* Explicit Scheme */
    /*
    The first and the last steps are separate.
    */
    for (auto m = 0; m < M; m++)
    {
        ump1[0] = u[0] + cfl * (u[1] - 2 * u[0] + u[J]);
        for (auto j = 1; j < J; j++)
        {
            ump1[j] = u[j] + cfl * (u[j + 1] - 2 * u[j] + u[j - 1]);
        }
        ump1[J] = u[J] + cfl * (u[0] - 2 * u[J] + u[J - 1]);
        plotVec.assign(ump1, ump1 + J + 1);
        matplot::plot(x, plotVec);
        matplot::xlim({XLIM_L, XLIM_R});
        matplot::ylim({YLIM_D, YLIM_U});

    }

    delete [] u;
    delete [] ump1;
    
    matplot::show();
}
