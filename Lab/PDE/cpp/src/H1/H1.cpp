#include "H1.hpp"
#include <vector>
#include <functional>
#include <iostream>
#include "../solver.hpp"


static void ode1EvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, double T, unsigned int M, double** &u, std::vector<double> &x);

int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        return -(x - 1) * (x + 1);
    };

    Solver odeSolver(u_initial);
    odeSolver.setScheme(ode1EvolveExplicitEulerPeriodic);
    odeSolver.solve();
    odeSolver.exportSolution("test.json");

    return 0;
}

static void ode1EvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, double T, unsigned int M, double** &u, std::vector<double> &x)
{
    double deltaT = (double) T / M;
    double deltaX = (b - a) / J;

    /* Initial Datum */
    for (auto j = 0; j <= J; j++)
    {
        x.push_back(a + (double) j * deltaX);
        u[0][j] = u_initial(x[j]);
    }

    /* Explicit Scheme */
    for (auto m = 0; m < M; m++)
    {
        for (auto j = 0; j <= J; j++)
        {
            u[m + 1][j] = u[m][j] * (1 - deltaT);
        }

        std::cout << "Progress: " << m << "/" << M << "    " << (double) m / M * 100 << std::endl;
    }

    
}
