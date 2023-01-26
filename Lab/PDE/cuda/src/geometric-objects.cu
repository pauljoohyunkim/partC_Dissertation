#include "geometric-objects.hpp"

/* Constructor for cuCurve */ 
cuCurve::cuCurve(unsigned int aJ)
{
    J = aJ;
    x.reserve(J);
    y.reserve(J);
    z.reserve(J);
    for (unsigned int i = 0; i < J; i++)
    {
        x.push_back(0.0);
        y.push_back(0.0);
        z.push_back(0.0);
    }
}

cuCurve::cuCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ)
{
    J = aX.size();
    x.reserve(J);
    y.reserve(J);
    z.reserve(J);
    for (unsigned int i = 0; i < J; i++)
    {
        x.push_back(aX[i]);
        y.push_back(aY[i]);
        z.push_back(aZ[i]);
    }
}
