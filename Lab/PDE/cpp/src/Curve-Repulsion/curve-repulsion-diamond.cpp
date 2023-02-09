#include "../math-objects.hpp"
#include "../geometric-objects.hpp"
#include "../solver.hpp"
#include "../export.hpp"
#include <vector>
#include <iostream>
#include <cmath>

#define LAMBDA 0.01
#define DELTAX 0.05
#define DELTAT 0.01
#define AZIMUTHAL_SPEED 0.5
#define ELEVATION 3
#define M 10000
#define PLOT_FREQUENCY 1

#define PI 3.14159265

static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, SolverCurveRepulsion& dk);

int main()
{
    SolverCurveRepulsion dk = SolverCurveRepulsion(2, 4, kernel1);

    /* Diamond */
    Vector3D x1{ 1.0, 0, 0 };
    Vector3D x2{ 0, 1.0, 0 };
    Vector3D x3{ -1.0, 0, 0 };
    Vector3D x4{ 0, -1.0, 0 };
    std::vector<Vector3D> veclist { x1, x2, x3, x4 };

    Curve c(veclist);
    auto J = c.getNPoints();

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

    /* Export to json */
    std::ofstream jsonX("x.json");
    std::ofstream jsonY("y.json");
    std::ofstream jsonZ("z.json");
    jsonX << "[";
    jsonY << "[";
    jsonZ << "[";
    


    /* Gradient Descent */
    for (auto t = 0; t < M; t++)
    {
        for (auto i = 0; i < J; i++)
        {
            d[i][0] = c[i][0] - dk.energyDifferential(c, Vector3D(DELTAX, 0, 0), i) / DELTAX * DELTAT - DELTAT * LAMBDA * c[i][0];
            d[i][1] = c[i][1] - dk.energyDifferential(c, Vector3D(0, DELTAX, 0), i) / DELTAX * DELTAT - DELTAT * LAMBDA * c[i][1];
            d[i][2] = c[i][2] - dk.energyDifferential(c, Vector3D(0, 0, DELTAX), i) / DELTAX * DELTAT - DELTAT * LAMBDA * c[i][2];
            x[i] = d[i][0];
            y[i] = d[i][1];
            z[i] = d[i][2];
        }
        x[J] = d[J][0];
        y[J] = d[J][1];
        z[J] = d[J][2];

        c = d;
        std::cout << "Progress: " << t << "/" << M << " (" << (float) t / M * 100 << "%)" << std::endl;
        if (t % PLOT_FREQUENCY == 0)
        {
            if (t != 0)
            {
                jsonX << ",\n";
                jsonY << ",\n";
                jsonZ << ",\n";
            }
            vectorParse(jsonX, x, J);
            vectorParse(jsonY, y, J);
            vectorParse(jsonZ, z, J);
        }
    }

    jsonX << "]";
    jsonY << "]";
    jsonZ << "]";
    jsonX.close();
    jsonY.close();
    jsonZ.close();

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
