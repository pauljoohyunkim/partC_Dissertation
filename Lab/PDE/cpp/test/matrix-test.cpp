#include "../src/math-objects.hpp"
#include <vector>
#include <iostream>

int main()
{

    Vector u(4);
    Vector v(4);

    u[0] = 1;
    u[1] = 2;
    u[2] = 3;
    v[1] = 10;

    auto w = u + v;



    return 0;
}
