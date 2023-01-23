#include "geometric-objects.hpp"

/* Curve Class */
/* Constructor
 * The parameter is the index of the last element.
 * Since 0 to M, there are M + 1 elements. */
Curve::Curve(unsigned int aM)
{
    points.reserve(M + 1);
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
