#include "../math-objects.hpp"
#include "../geometric-objects.hpp"
#include "discrete-kernel.hpp"
#include "curve-repulsion-1.hpp"
#include <vector>
#include <iostream>
#include <matplot/matplot.h>

#define LAMBDA 0.01

static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, DiscreteKernel& dk);

int main()
{
    DiscreteKernel dk = DiscreteKernel(2, 4, kernel1);

    /* Points on the curve (or rather... polygon)
     * x1(2.5, 3, 2)
     * x2(-3.5, 3.5, -2)
     * x3(-3.2, -3.25, 1.8)
     * x4(5.5, -6.5, 0)
     * */

    Vector3D x1(13.2, -16.4, 11.12);
    Vector3D x2(10.18, 9.09, 0);
    Vector3D x3(16.53, -25.6, -6.01);
    Vector3D x4(-5.91, -21.32, 0);

    std::vector<Vector3D> veclist { x1, x2, x3, x4 };

    Curve c(veclist);
    auto J = c.getNPoints();

    double energy = dk.energy(c);
    double energydiff = dk.energyDifferential(c, Vector3D(0.1, 0, 0), 1);

    Curve d(c.getNPoints());


    /* For plotting */
    std::vector<double> x {} ;
    std::vector<double> y {} ;
    std::vector<double> z {} ;
    x.reserve(J + 1);
    y.reserve(J + 1);
    z.reserve(J + 1);
    for (auto i = 0; i < J + 1; i++)
    {
        x.push_back(0);
        y.push_back(0);
        z.push_back(0);
    }

    


    /* Gradient Descent */
    for (auto t = 0; t < 400; t++)
    {
        for (auto i = 0; i < J; i++)
        {
            d[i][0] = c[i][0] - dk.energyDifferential(c, Vector3D(0.1, 0, 0), i) - LAMBDA * c[i][0];
            d[i][1] = c[i][1] - dk.energyDifferential(c, Vector3D(0, 0.1, 0), i) - LAMBDA * c[i][1];
            d[i][2] = c[i][2] - dk.energyDifferential(c, Vector3D(0, 0, 0.1), i) - LAMBDA * c[i][2];
            x[i] = d[i][0];
            y[i] = d[i][1];
            z[i] = d[i][2];
        }
        x[J] = d[J][0];
        y[J] = d[J][1];
        z[J] = d[J][2];

        c = d;
        std::cout << t << ": " << dk.energy(c) << std::endl;
        matplot::plot3(x, y, z);
        //matplot::xrange({-30, 30});
        //matplot::yrange({-30, 30});

        matplot::save(std::to_string(t) + ".png");
    }

    return 0;
}

/* Based on "Trapezoid rule":
 * Takes average based on four points
 * */
static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, DiscreteKernel& dk)
{
    double kij { 0 };  

    kij += dk.kernelalphabeta(xi, xj, Ti);
    kij += dk.kernelalphabeta(xi, xjp1, Ti);
    kij += dk.kernelalphabeta(xip1, xj, Ti);
    kij += dk.kernelalphabeta(xip1, xjp1, Ti);

    return kij / 4;

}
