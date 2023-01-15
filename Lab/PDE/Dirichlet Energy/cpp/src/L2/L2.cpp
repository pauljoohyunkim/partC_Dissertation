#include <iostream>
#include <functional>
#include <algorithm>
#include <iterator>
#include <vector>
#include <matplot/matplot.h>
#include "../solver.hpp"
#include "../params.hpp"
#include "L2.hpp"

static void heatEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, unsigned int T, unsigned int M, double** &u);

int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        return -(x - 1) * (x + 1);
    };

    Solver heatEquationSolver(u_initial);
    heatEquationSolver.setScheme(heatEvolveExplicitEulerPeriodic);
    heatEquationSolver.solve();

    return 0;
}


static void heatEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, unsigned int T, unsigned int M, double** &u)
{
    double deltaT = (double) T / M;
    double deltaX = (b - a) / J;
    double cfl = deltaT / (deltaX * deltaX);
    std::vector<double> x {};
    
    if (cfl > 0.5)
    {
        std::cout << "Warning: the CFL number is " << cfl << ", so stability is not guaranteed";
    }

    /* Initial Datum */
    for (auto j = 0; j <= J; j++)
    {
        x.push_back(a + (double) j * deltaX);
        u[0][j] = u_initial(x[j]);
    }


    /* Explicit Scheme */
    /*
    The first and the last steps are separate.
    */
    for (auto m = 0; m < M; m++)
    {
        u[m + 1][0] = u[m][0] + cfl * (u[m][1] - 2 * u[m][0] + u[m][J]);
        for (auto j = 1; j < J; j++)
        {
            u[m + 1][j] = u[m][j] + cfl * (u[m][j + 1] - 2 * u[m][j] + u[m][j - 1]);
        }
        u[m + 1][J] = u[m][J] + cfl * (u[m][0] - 2 * u[m][J] + u[m][J - 1]);

        std::cout << "Progress: " << m << "/" << M << "    " << (double) m / M * 100 << std::endl;
    }
}
