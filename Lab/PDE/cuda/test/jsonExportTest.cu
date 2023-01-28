#include "../src/export.hpp"
#include <iostream>
#include <vector>

int main()
{
    std::vector<double> v { 1, 2, 3, 4 };

    vectorParse(std::cout, v, 4);

    return 0;
}
