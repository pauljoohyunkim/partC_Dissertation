#include "math-objects.hpp"

/* Constructor for Matrix */
Matrix::Matrix(unsigned int m, unsigned int n)
{
    dimensions.first = m;
    dimensions.second = n;
}

unsigned int Matrix::getNRows()
{
    return dimensions.first;
}

unsigned int Matrix::getNColumns()
{
    return dimensions.second;
}

std::pair<unsigned int, unsigned int> Matrix::getNRowsColumns()
{
    return dimensions;
}
