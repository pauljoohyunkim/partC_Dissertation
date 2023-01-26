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

    std::cout << C.getValFromDevice(0, -1) << std::endl;
    std::cout << C.getValFromDevice(1, -1) << std::endl;
    std::cout << C.getValFromDevice(2, -1) << std::endl;
    std::cout << C.getValFromDevice(0, 1) << std::endl;
    std::cout << C.getValFromDevice(1, 1) << std::endl;
    std::cout << C.getValFromDevice(2, 1) << std::endl;

    return 0;
}
