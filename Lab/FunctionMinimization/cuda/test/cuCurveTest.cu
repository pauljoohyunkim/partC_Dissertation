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

    curve.cudafy();

    //printCoefficientsPartiallyDEBUG<<<1,1>>>(curve.dev_cos_table + 0);
    //printCoefficientsPartiallyDEBUG<<<1,1>>>(curve.dev_cos_table + 1);
    //printCoefficientsPartiallyDEBUG<<<1,1>>>(curve.dev_cos_table + 2);
    //printCoefficientsPartiallyDEBUG<<<1,1>>>(curve.dev_cos_table + 3);
    //printCoefficientsPartiallyDEBUG<<<1,1>>>(curve.dev_cos_table + 4);
    //printCoefficientsPartiallyDEBUG<<<1,1>>>(curve.dev_cos_table + 5);
    //printCoefficientsPartiallyDEBUG<<<1,1>>>(&curve.dev_cos_table[2 + 3 * 1]);
    //queryDEBUG<<<1,1>>>(curve.dev_cos_table, 1, 2, curve.J);

    fillDEBUG<<<1,1>>>(curve.dev_x, curve.dev_y, curve.dev_z, curve.dev_coefficients, curve.dev_cos_table, curve.dev_sin_table, curve.resolution, curve.J);

    dim3 grid(curve.resolution, curve.resolution);
    tangentPointEnergyMatrixFill<<<grid, 1>>>(curve.dev_x, curve.dev_y, curve.dev_z, curve.dev_energy_matrix, curve.resolution);

    curve.cudaFlush();

    return 0;
}
