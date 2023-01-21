#include "../src/math-objects.hpp"
#include <vector>
#include <iostream>

int main()
{
    double** a = new double* [2];
    a[0] = new double [2];
    a[1] = new double [2];


    double** b = new double* [2];
    b[0] = new double [2];
    b[1] = new double [2];

    a[0][0] = -2;
    a[0][1] = 1;
    a[1][0] = 0;
    a[1][1] = 4;

    b[0][0] = 6;
    b[0][1] = 5;
    b[1][0] = -7;
    b[1][1] = 1;

    Matrix A(2, 2, a);
    Matrix B(2, 2, b);

    auto C = A + B;

    std::cout << C[0][1] << std::endl;

    auto v4 = std::vector<double> {1,2,3};
    Vector3D v(v4);

    return 0;
}
