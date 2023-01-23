#ifndef GEOMETRIC_OBJECTS_HPP
#define GEOMETRIC_OBJECTS_HPP

#include "math-objects.hpp"
#include <vector>

class Curve
{
    public:
        /* Constructor */
        Curve(unsigned int M);


        /* Get i-th point. Loop around the index.
         * -1 for the last
         *  M for the first */
        operator [] (int i);
    protected:
        std::vector<Vector3D> points;
}

#endif  // geometric-objects.hpp
