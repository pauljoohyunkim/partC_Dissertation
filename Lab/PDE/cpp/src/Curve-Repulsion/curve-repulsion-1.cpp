#include "../math-objects.hpp"
#include "discrete-kernel.hpp"

static double kernel1(Vector3D& xi, Vector3D& xip1, Vector3D& xj, Vector3D& xjp1, Vector3D& Ti, DiscreteKernel& dk);

int main()
{
    DiscreteKernel dk = DiscreteKernel(2, 4, kernel1);
    Vector3D v1(0, 0, 1);
    Vector3D v2(0, 1, 0);
    Vector3D v3(1, 0, 2);
    Vector3D v4(-1, 2, 1);
    Vector3D v5(-3, 2, 1);

    auto val = kernel1(v1, v2, v3, v4, v5, dk);

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
