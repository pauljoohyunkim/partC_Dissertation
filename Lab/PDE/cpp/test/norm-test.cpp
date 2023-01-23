#include "../src/math-objects.hpp"
#include <iostream>

int main()
{
    Vector3D v(3, 4, 5);

    std::cout << l2norm(v) << std::endl;

    return 0;

}
