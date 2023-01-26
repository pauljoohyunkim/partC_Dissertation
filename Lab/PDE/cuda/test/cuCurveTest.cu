#include "../src/geometric-objects.hpp"
#include <vector>

int main()
{
    std::vector<double> x = { 1, 2, 3, 4, 5, 6 };
    std::vector<double> y = { 0, 2, 4, 6, 8, -1 };
    std::vector<double> z = { -1, -2, 3, -4, 5, 0.6 };
    cuCurve C(x, y, z);
    C.cudafy();

    return 0;
}
