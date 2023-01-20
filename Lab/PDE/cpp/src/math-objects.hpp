#ifndef MATH_OBJECTS_HPP
#define MATH_OBJECTS_HPP

#include <utility>

class Matrix
{
    public:
        /* Matrix Constructor */
        Matrix(unsigned int m, unsigned int n);

        unsigned int getNRows();
        unsigned int getNColumns();
        std::pair<unsigned int, unsigned int> getNRowsColumns();

    private:
        std::pair<unsigned int, unsigned int> dimensions {1, 1};


};


#endif  // MATH_OBJECTS_HPP
