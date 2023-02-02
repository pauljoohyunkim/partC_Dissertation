#include "../src/solver.hpp"
#include "../src/geometric-objects.hpp"
#include "../src/math-objects.hpp"
#include <vector>
#include <iomanip>


static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, SolverCurveRepulsion& dk);

int main()
{
    SolverCurveRepulsion dk = SolverCurveRepulsion(2, 4, kernel1);
    Vector3D x1{ 1, 0, 0 };
    Vector3D x2{ 0, 1, 0 };
    Vector3D x3{ -1, 0, 0 };
    Vector3D x4{ 0, -1, 0 };
    std::vector<Vector3D> veclist { x1, x2, x3, x4 };

    Curve c(veclist);
    auto J = c.getNPoints();

    auto enrgy = dk.energy(c);

    std::cout << std::setprecision(15) << enrgy << std::endl;
    
}

static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, SolverCurveRepulsion& dk)
{
    double kij { 0 };  

    kij += dk.kernelalphabeta(xi, xj, Ti);
    kij += dk.kernelalphabeta(xi, xjp1, Ti);
    kij += dk.kernelalphabeta(xip1, xj, Ti);
    kij += dk.kernelalphabeta(xip1, xjp1, Ti);

    return kij / 4;

}
