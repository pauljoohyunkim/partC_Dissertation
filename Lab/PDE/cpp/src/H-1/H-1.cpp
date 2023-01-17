#include "H-1.hpp"
#include "../solver.hpp"
#include <vector>
#include <iostream>
#include <functional>

static void biharmonicEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, unsigned int T, unsigned int M, double** &u, std::vector<double> &x);


int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        return -10 * (x - 1) * (x + 1);
    };

    Solver biharmonicEquationSolver(u_initial, -1, 1, 20, 10, 1000000);
    biharmonicEquationSolver.setScheme(biharmonicEvolveExplicitEulerPeriodic);
    biharmonicEquationSolver.solve();
    biharmonicEquationSolver.exportSolution("test.json", 10);
    
    return 0;
}

static void biharmonicEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a, double b, unsigned int J, unsigned int T, unsigned int M, double** &u, std::vector<double> &x)
{
    double deltaT = (double) T / M;
    double deltaX = (b - a) / J;
    double cfl = deltaT / (deltaX * deltaX * deltaX * deltaX);

    if (cfl > 0.125)
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
    The first two and the last two steps are separate.
    */
    for (auto m = 0; m < M; m++)
    {
        u[m + 1][0] = u[m][0] - cfl * (u[m][2] - 4 * u[m][1] + 6 * u[m][0] - 4 * u[m][J] + u[m][J - 1]);
        u[m + 1][1] = u[m][1] - cfl * (u[m][3] - 4 * u[m][2] + 6 * u[m][1] - 4 * u[m][0] + u[m][J]);
        for (auto j = 2; j < J - 1; j++)
        {
            u[m + 1][j] = u[m][j] - cfl * (u[m][j + 2] - 4 * u[m][j + 1] + 6 * u[m][j] - 4 * u[m][j - 2] + u[m][j - 2]);
        }
        u[m + 1][J - 1] = u[m][J - 1] - cfl * (u[m][0] - 4 * u[m][J] + 6 * u[m][J - 1] - 4 * u[m][J - 3] + u[m][J - 3]);
        u[m + 1][J] = u[m][J] - cfl * (u[m][1] - 4 * u[m][0] + 6 * u[m][J] - 4 * u[m][J - 1] + u[m][J - 3]);

        std::cout << "Progress: " << m << "/" << M << "    " << (double) m / M * 100 << std::endl;
    }

}
