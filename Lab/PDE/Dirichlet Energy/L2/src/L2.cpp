#include <iostream>
#include <matplot/matplot.h>
#include "solvers.hpp"
#include "L2.hpp"

int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        return -10 * (x - 1) * (x + 1);
    };

    heatEvolveExplicitEulerPeriodic(u_initial);
    
}
