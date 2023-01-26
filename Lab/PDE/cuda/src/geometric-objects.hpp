#ifndef GEOMETRIC_OBJECTS_H
#define GEOMETRIC_OBJECTS_H

#include <vector>

struct cuCurve
{
    public:
        /* Constructor */
        cuCurve(unsigned int aJ);
        cuCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ);

        /* Variables */
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        unsigned int J;

};


#endif  // geometric-objects.h
