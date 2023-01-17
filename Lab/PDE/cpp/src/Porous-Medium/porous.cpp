#include "porous.hpp"
#include "../solver.hpp"
#include <iostream>
#include <functional>
#include <vector>


static void porousMediumExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, double T, unsigned int M, double** &u, std::vector<double> &x);

int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        if(x < -0.2 || x > 0.2)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    };

    Solver porousMediumEquationSolver(u_initial, -1, 1, 400, 0.25, 80000);
    porousMediumEquationSolver.setScheme(porousMediumExplicitEulerPeriodic);
    porousMediumEquationSolver.solve();
    porousMediumEquationSolver.exportSolution("test.json", 10);

    return 0;
}

static void porousMediumExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, double T, unsigned int M, double** &u, std::vector<double> &x)
{
    double deltaT = (double) T / M;
    double deltaX = (b - a) / J;
    double cfl = deltaT / (deltaX * deltaX);

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
        u[m + 1][0] = u[m][0] + cfl * (u[m][1] * u[m][1] - 2 * (u[m][0] * u[m][0]) + (u[m][J] * u[m][J]));
        for (auto j = 1; j < J; j++)
        {
            u[m + 1][j] = u[m][j] + cfl * (u[m][j + 1] * u[m][j + 1] - 2 * (u[m][j] * u[m][j]) + (u[m][j - 1] * u[m][j - 1]));
        }
        u[m + 1][J] = u[m][J] + cfl * (u[m][0] * u[m][0] - 2 * (u[m][J] * u[m][J]) + (u[m][J - 1] * u[m][J - 1]));

        std::cout << "Progress: " << m << "/" << M << "    " << (double) m / M * 100 << std::endl;
    }


}

