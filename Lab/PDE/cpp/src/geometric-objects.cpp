#include "geometric-objects.hpp"

/* Curve Class */
/* Constructor
 * The parameter is the number of points. */
Curve::Curve(unsigned int aM)
{
    points.reserve(M);
    M = aM;
}

Curve::Curve(std::vector<Vector3D> &veclist)
{
    points = veclist;
    M = points.size();
}

Vector3D& Curve::operator [] (int i)
{
    return points[((i % (int) M) + M) % M];
}
