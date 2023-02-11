#include "../src/curve.hpp"
#include "../src/curverepulsion.hpp"
#include <vector>

int main()
{
    std::vector<double> xa { 1, 2, 3 };
    std::vector<double> xb { 4, 5, 6 };
    std::vector<double> ya { 7, 8, 9 };
    std::vector<double> yb { 10, 11, 12 };
    std::vector<double> za { 13, 14, 15 };
    std::vector<double> zb { 16, 17, 18 };

    FourierCurve curve(xa, xb, ya, yb, za, zb);


    fillDifferentialMatrix(curve, 0.1);

    gradientDescent<<<6 * (curve.J + 1), 1>>>(curve.dev_coefficients, curve.dev_differential_coefficients, 0.1, curve.J);

    curve.cudaFlush();


    return 0;
}
