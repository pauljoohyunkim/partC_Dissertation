#include "../math-objects.hpp"
#include "../geometric-objects.hpp"
#include "discrete-kernel.hpp"
#include "curve-repulsion-1.hpp"
#include <vector>

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

    Vector3D x1(2.5, 3, 2);
    Vector3D x2(-3.5, 3.5, -2);
    Vector3D x3(-3.2, -3.25, 1.8);
    Vector3D x4(5.5, -6.5, 0);

    std::vector<Vector3D> veclist { x1, x2, x3, x4 };

    Curve c(veclist);

    double energy = dk.energy(c);

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
