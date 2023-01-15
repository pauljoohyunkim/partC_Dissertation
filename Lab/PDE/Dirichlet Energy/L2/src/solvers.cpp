#include "solvers.hpp"
#include <functional>
#include <iostream>
#include <vector>
#include <matplot/matplot.h>

void heatEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, unsigned int T, unsigned int M)
{
    double deltaT = (double) T / M;
    double deltaX = (b - a) / J;
    double cfl = deltaT / (deltaX * deltaX);
    std::vector<double> u {};
    std::vector<double> x {};                   /* Vector of x values */
    std::vector<double> ump1 {};
    
    if (cfl > 0.5)
    {
        std::cout << "Warning: the CFL number is " << cfl << ", so stability is not guaranteed";
    }

    /* Initial Datum */
    for (auto j = 0; j <= J; j++)
    {
        x.push_back(a + (double) j * deltaX);
        u.push_back(u_initial(x.back()));
    }


    matplot::hold(matplot::on);
    
    /* Explicit Scheme */
    /*
    The first and the last steps are separate.
    */
    for (auto m = 0; m < M; m++)
    {
        ump1.clear();
        ump1.push_back(u[0] + cfl * (u[1] - 2 * u[0] + u[J]));
        for (auto j = 1; j < J; j++)
        {
            ump1.push_back(u[j] + cfl * (u[j + 1] - 2 * u[j] + u[j - 1]));
        }
        ump1.push_back(u[J] + cfl * (u[0] - 2 * u[J] + u[J - 1]));
        u = ump1;
        matplot::plot(x, ump1);
    }
    
    matplot::show();
}
