#include "../src/math-objects.hpp"
#include <vector>
#include <iostream>

int main()
{
    double** a = new double* [3];
    a[0] = new double [3];
    a[1] = new double [3];
    a[2] = new double [3];


    double** b = new double* [3];
    b[0] = new double [3];
    b[1] = new double [3];
    b[2] = new double [3];

    a[0][0] = -2;
    a[0][1] = 1;
    a[1][0] = 0;
    a[1][1] = 4;

    b[0][0] = 6;
    b[0][1] = 5;
    b[1][0] = -7;
    b[1][1] = 1;

    Matrix A(3, 3, a);
    //Matrix B(3, 3, b);

    /*
    auto C = A + B;

    std::cout << C[0][1] << std::endl;
    */

    Vector3D v1(2,-3,-1);
    Vector3D v2(1,4,-2);

    /*
    auto v3 = v1 ^ v2;
    auto x = v1 % v2;
    auto vA = vectorize(A);
    auto vA3 = vectorize(vA);
    */

    auto B = A * v1;
    auto C = A * v2;

    auto Av1 = vectorize(B);
    auto Av2 = vectorize3d(C);

    return 0;
}
