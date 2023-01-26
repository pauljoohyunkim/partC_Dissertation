#include "solver.hpp"
//#include "geometric-objects.hpp"

__device__ void cross (double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3)
{
    x3 = y1 * z2 - y2 * z1;
    y3 = z1 * x2 - x1 * z2;
    z3 = x1 * y2 - x2 * y1;
}
