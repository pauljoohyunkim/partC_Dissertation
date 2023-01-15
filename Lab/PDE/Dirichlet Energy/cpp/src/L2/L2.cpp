#include <iostream>
#include <matplot/matplot.h>
#include "../solver.hpp"
#include "L2.hpp"

int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        return -(x - 1) * (x + 1);
    };

    Solver heatEquationSolver(u_initial);

    return 0;
}
