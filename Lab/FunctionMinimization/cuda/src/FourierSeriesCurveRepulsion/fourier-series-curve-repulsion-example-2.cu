#include "../curve.hpp"
#include "../curverepulsion.hpp"
#include "../export.hpp"
#include <vector>

#define CURVE_RESOLUTION 30
#define PERTURBATION 0.01
#define STEPSIZE 0.01
#define M 10000

int main()
{
    /* Predefined curve */
    std::vector<double> xa { -0.63538135, -0.30286386, -0.54336102,  1.90359802,  1.78795006,
        1.79616879,  0.9555309 , -0.28603216,  1.87019586,  0.25694016};
    std::vector<double> xb { 0.51252115, -1.43472404, -1.51223027, -1.1053243 ,  1.68650192,
        -0.54617886, -0.91839065, -0.5668194 , -0.69252904,  1.82442381 };
    std::vector<double> ya {-0.81930084, -0.75651476,  0.87815106,  0.51176656,  1.50663107,
        -1.20918657,  0.52859345, -0.46538907, -1.71020548,  0.94714671 };
    std::vector<double> yb {0.28624529, -1.51078498,  1.76745975, -0.07898508,  1.492646  ,
        -1.33914639, -0.87156463, -0.15177939,  1.17117152,  0.39036983};
    std::vector<double> za {0.46658214,  0.25219143, -1.22422532, -1.03771514, -0.46819253,
        -0.20119577, -1.37406833, -0.33057391,  0.4935988 , -1.31441903};
    std::vector<double> zb {0.82400325,  0.14854623, -1.52653151,  0.11398457,  0.35415078,
        -0.03107083, -0.83402719, -0.3050666 , -1.01450412, -1.75535255};
    FourierCurve curve(xa, xb, ya, yb, za, zb, CURVE_RESOLUTION);


    JsonExporter xExporter { "x.json" };
    JsonExporter yExporter { "y.json" };
    JsonExporter zExporter { "z.json" };
    

    for (unsigned int t = 0; t < M; t++)
    {
        fillDifferentialMatrix(curve, PERTURBATION);
        gradientDescent<<<6 * (curve.J + 1), 1>>>(curve.dev_coefficients, curve.dev_differential_coefficients, STEPSIZE, curve.J);
        curve.cudaFlush();

        xExporter << curve.x;
        yExporter << curve.y;
        zExporter << curve.z;

        std::cout << "Progress: " << t << "/" << M << std::endl;
    }

    return 0;
}
