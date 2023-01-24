#include "geometric-objects.hpp"

/* Curve Class */
/* Constructor
 * The parameter is the number of points. */
Curve::Curve(unsigned int aJ)
{
    points.reserve(J);
    J = aJ;
}

Curve::Curve(std::vector<Vector3D> &veclist)
{
    points = veclist;
    J = points.size();
}

Vector3D& Curve::operator [] (int i)
{
    return points[((i % (int) J) + J) % J];
}

unsigned int Curve::getNPoints()
{
    return J;
}
