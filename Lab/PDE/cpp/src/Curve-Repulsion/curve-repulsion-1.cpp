#include "../math-objects.hpp"
#include "../geometric-objects.hpp"
#include "../solver.hpp"
#include "curve-repulsion-1.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <matplot/matplot.h>

#define LAMBDA 0.01
#define DELTAX 0.01
#define DELTAT 0.005
#define AZIMUTHAL_SPEED 0.5
#define ELEVATION 3
#define M 10000
#define PLOT_FREQUENCY 2

#define PI 3.14159265

static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, SolverCurveRepulsion& dk);

int main()
{
    SolverCurveRepulsion dk = SolverCurveRepulsion(2, 4, kernel1);

    /* Points on the curve (or rather... polygon)
     * x1(2.5, 3, 2)
     * x2(-3.5, 3.5, -2)
     * x3(-3.2, -3.25, 1.8)
     * x4(5.5, -6.5, 0)
     * */

    /* Example 1 */
    //Vector3D x1(13.2, -16.4, 11.12);
    //Vector3D x2(10.18, 9.09, 0);
    //Vector3D x3(16.53, -25.6, -6.01);
    //Vector3D x4(-5.91, -21.32, 0);
    //std::vector<Vector3D> veclist { x1, x2, x3, x4 };

    /* Example 2 */
    //Vector3D x1{ -2.75, 2.97, 0 };
    //Vector3D x2{ -0.6, -2.38, 2 };
    //Vector3D x3{ 3.38, -3.62, -1.55 };
    //Vector3D x4{ 3.1, 2.23, 0 };
    //Vector3D x5{ 0.74, 3.36, 3.24 };
    //Vector3D x6{ -6.71, 1.73, 4 };
    //Vector3D x7{ 1.75, -6.61, -2.28 };
    //Vector3D x8{ 2.6, 3.21, 5.24 };
    //Vector3D x9{ -4.7, 6.23, 4.38 };
    //std::vector<Vector3D> veclist { x1, x2, x3, x4, x5, x6, x7, x8, x9 };
    //

    Vector3D x1 {1,0,-1};
    Vector3D x2 {2,2,-2};
    Vector3D x3 {3,4,3};
    Vector3D x4 {4,6,-4};
    Vector3D x5 {5,8,5};
    Vector3D x6 {6,-1,0.6};
    /* Example 3: Helix + Semicircle */
    //const int resolution = 16;
    std::vector<Vector3D> veclist { x1, x2, x3, x4, x5, x6 };
    //for (auto i = 0; i < resolution; i++)
    //{
    //    double theta = 4 * PI * (double) i / resolution;
    //    Vector3D p(cos(theta), sin(theta), theta / (2 * PI));
    //    veclist.push_back(p);
    //}
    //for (auto i = 1; i < resolution; i++)
    //{
    //    double theta = PI * (double) i / resolution;
    //    Vector3D p(1, 2 * sin(theta), 1 + cos(theta));
    //    veclist.push_back(p);
    //}

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
    for (auto t = 0; t < M; t++)
    {
        for (auto i = 0; i < J; i++)
        {
            d[i][0] = c[i][0] - dk.energyDifferential(c, Vector3D(DELTAX, 0, 0), i) / DELTAX * DELTAT - DELTAT * LAMBDA * (2 * c[i][0] - c[i + 1][0] - c[i - 1][0]);
            d[i][1] = c[i][1] - dk.energyDifferential(c, Vector3D(0, DELTAX, 0), i) / DELTAX * DELTAT - DELTAT * LAMBDA * (2 * c[i][1] - c[i + 1][1] - c[i - 1][1]);
            d[i][2] = c[i][2] - dk.energyDifferential(c, Vector3D(0, 0, DELTAX), i) / DELTAX * DELTAT - DELTAT * LAMBDA * (2 * c[i][2] - c[i + 1][2] - c[i - 1][2]);
            x[i] = d[i][0];
            y[i] = d[i][1];
            z[i] = d[i][2];
        }
        x[J] = d[J][0];
        y[J] = d[J][1];
        z[J] = d[J][2];

        c = d;
        std::cout << t << ": " << dk.energy(c) << std::endl;
        if (t % PLOT_FREQUENCY == 0)
        {
            auto curvePlot = matplot::plot3(x, y, z);
            curvePlot->line_width(5);
            matplot::view(AZIMUTHAL_SPEED * t, ELEVATION);
            matplot::xrange({-5, 5});
            matplot::yrange({-5, 5});
            //matplot::show();

            //matplot::save(std::to_string(t) + ".png");
        }
        //if (t == M-2)
        //{
        //    matplot::show();
        //}
    }

    return 0;
}

/* Based on "Trapezoid rule":
 * Takes average based on four points
 * */
static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, SolverCurveRepulsion& dk)
{
    double kij { 0 };  

    kij += dk.kernelalphabeta(xi, xj, Ti);
    kij += dk.kernelalphabeta(xi, xjp1, Ti);
    kij += dk.kernelalphabeta(xip1, xj, Ti);
    kij += dk.kernelalphabeta(xip1, xjp1, Ti);

    return kij / 4;

}
