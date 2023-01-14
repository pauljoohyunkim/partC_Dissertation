#include "solvers.hpp"
#include <functional>
#include <iostream>
#include <vector>
#include <matplot/matplot.h>

void heatEvolveExplicitEulerPeriodic(double a = -1, double b = 1, unsigned int J = 20, unsigned int T = 10000, unsigned int M = 10000000, std::function<double(double)> u_initial)
{
    auto deltaT = T / M;
    auto deltaX = (b - a) / J;
    auto cfl = deltaT / (deltaX * deltaX);
    std::vector<double> u {};
    std::vector<double> ump1 {};
    
    if (cfl > 0.5)
    {
        std::cout << "Warning: the CFL number is " << cfl << ", so stability is not guaranteed";
    }

    /* Initial Datum */
    for (auto j = 0; j < J; j++)
    {
        u.push_back(u_initial(a + (double) j * deltaX));
    }
    
    /* Explicit Scheme */
    /*
    The first and the last steps are separate.
    */
    for (auto m = 0; m < M, m++)
    {
        ump1.clear();
        ump1.push_back(u[0] + cfl * (u[1] - 2 * u[0] + u[J]));
        for (auto j = 1; j < J; j++)
        {
            ump1.push_back(u[j] + cfl * (u[j + 1] - 2 * u[j] + u[j - 1]));
        }
        ump1.push_back(u[J] + cfl * (u[0] - 2 * u[J] + u[J - 1]));
        u = ump1;
    }


}
