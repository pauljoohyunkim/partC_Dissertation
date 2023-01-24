#ifndef GEOMETRIC_OBJECTS_HPP
#define GEOMETRIC_OBJECTS_HPP

#include "math-objects.hpp"
#include <vector>

class Curve
{
    public:
        /* Constructor */
        Curve(unsigned int aJ);
        Curve(std::vector<Vector3D> &veclist);


        /* Get i-th point. Loop around the index.
         * -1 for the last
         *  M for the first */
        Vector3D& operator [] (int i);
        /* Get number of points */
        unsigned int getNPoints();
    protected:
        std::vector<Vector3D> points;
        unsigned int J { 0 };
};

#endif  // geometric-objects.hpp
