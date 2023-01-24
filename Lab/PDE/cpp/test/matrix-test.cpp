#include "../src/math-objects.hpp"
#include <vector>
#include <iostream>

int main()
{

    double** b = new double* [3];
    for (auto i = 0; i < 3; i++)
    {
        b[i] = new double [3];
    }

    b[0][0] = -2;
    b[0][1] = 2;
    b[1][1] = 1;
    b[2][1] = 3;

    Matrix B(3, 3, b);

    Vector3D v(1);

    auto Bv = matvecmul3D(B, v);
    auto Bv3 = Bv * 3;

    auto B3 = B * 3.0;

    for (auto i = 0; i < 3; i++)
    {
        delete [] b[i];
    }
    delete [] b;

    return 0;
}
