#include "math-objects.hpp"

/* Constructor for Matrix */
Matrix::Matrix(unsigned int m, unsigned int n, double initval)
{
    dimensions.first = m;
    dimensions.second = n;
    rawMatrix.resize(m);
    for (unsigned int mi = 0; mi < m; mi++)
    {
        rawMatrix[mi].resize(n, initval);
    }
}

Matrix::Matrix(unsigned int m, unsigned int n, double** &pMatrix)
{
    dimensions.first = m;
    dimensions.second = n;
    rawMatrix.resize(m);
    /* Copy values from pMatrix */
    for (unsigned int mi = 0; mi < m; mi++)
    {
        for (unsigned int ni = 0; ni < n; ni++)
        {
            rawMatrix[mi].push_back(pMatrix[mi][ni]);
        }
    }

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
