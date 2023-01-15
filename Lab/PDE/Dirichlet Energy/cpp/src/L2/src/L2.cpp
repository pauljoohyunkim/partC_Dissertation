#include <iostream>
#include <matplot/matplot.h>
#include "solvers.hpp"
#include "L2.hpp"

int main()
{
    /* The initial datum */
    auto u_initial = [] (double x) -> double
    {
        return -(x - 1) * (x + 1);
    };

    heatEvolveExplicitEulerPeriodic(u_initial, -1, 1, 40, 10, 10000);
    
}
