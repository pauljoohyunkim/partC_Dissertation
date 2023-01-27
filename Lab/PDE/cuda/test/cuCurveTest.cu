#include "../src/geometric-objects.hpp"
#include "../src/solver.hpp"
#include <vector>
#include <iostream>

int main()
{
    std::vector<double> x = { 1, 2, 3, 4, 5, 6 };
    std::vector<double> y = { 0, 2, 4, 6, 8, -1 };
    std::vector<double> z = { -1, -2, 3, -4, 5, 0.6 };
    cuRepulsiveCurve C(x, y, z);
    C.cudafy();

    dim3 grid(C.J, C.J);
    fillEnergyMatrix<<<grid, 1>>>(C.dev_x, C.dev_y, C.dev_z, C.dev_energyMatrix, C.J);
    C.flushFromDevice();

    return 0;
}
