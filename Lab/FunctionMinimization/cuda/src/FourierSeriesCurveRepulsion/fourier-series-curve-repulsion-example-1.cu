#include "../curve.hpp"
#include "../curverepulsion.hpp"
#include "../export.hpp"
#include <vector>

#define CURVE_RESOLUTION 30
#define PERTURBATION 0.01
#define STEPSIZE 0.1
#define M 10000

int main()
{
    /* Predefined curve */
    std::vector<double> xa { -0.15272798, -1.60413542,  0.9939416 ,  2.95901542,  2.44886558,
                             -2.77151784, -1.66386979,  1.10095863,  1.09430543,  0.85264909,
                              0.18104473, -0.98786042,  0.66693053,  0.01043734,  0.01420098 };
    std::vector<double> xb { -1.35163415, -0.98074608, -2.14839559, -1.23042562,  2.15706847,
                              1.45403898, -0.00448405,  0.40229761,  0.61911565,  1.60795054,
                              1.47789099, -2.05937418,  1.09084598, -1.81701127, -1.4206581  };
    std::vector<double> ya {  0.57127358, -2.34858219,  0.903887  , -1.29340586, -2.40634033,
                              2.99691204, -2.40940774, -0.37178444,  0.32742878, -0.5715192 ,
                              1.1199972 ,  1.75841556,  0.68520055, -1.07542582, -1.73604004 };
    std::vector<double> yb { -2.50469606,  2.35970779, -0.88397131,  0.62306188, -2.57317773,
                              1.65548942, -2.35097927, -1.15938552, -2.951843  ,  2.42684191,
                             -0.95732087, -2.11943218, -1.86355292, -2.30007547,  0.91228002 };
    std::vector<double> za { -0.75211074, -2.03187897,  1.29089092, -0.11023411,  2.35830628,
                              1.665579  , -2.68856911, -0.16545203, -2.81192059,  1.58216797,
                              0.37327507, -0.12133424, -2.52793162,  0.11757547,  1.22681207 };
    std::vector<double> zb {  -0.93119624,  0.29444041,  0.4587438 , -1.63970205, -0.5378548 ,
                              0.29823082, -1.60149189,  2.40852357,  2.97347022, -1.80613509,
                              0.25954272,  2.47627695, -0.33145004, -1.0766391 , -1.52592399 };
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
