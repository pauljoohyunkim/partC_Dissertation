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
        rawMatrix[mi].resize(n);
        for (unsigned int ni = 0; ni < n; ni++)
        {
            rawMatrix[mi][ni] = pMatrix[mi][ni];
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

Matrix Matrix::operator + (Matrix &M)
{
    /* Follow the format from the current matrix. */
    Matrix matrix(this->dimensions.first, this->dimensions.second);
    for (unsigned int mi = 0; mi < matrix.getNRows(); mi++)
    {
        for (unsigned int ni = 0; ni < matrix.getNColumns(); ni++)
        {
            matrix.rawMatrix[mi][ni] = this->rawMatrix[mi][ni] + M.rawMatrix[mi][ni];
        }
    }

    return matrix;
}

Matrix Matrix::operator - (Matrix &M)
{
    /* Follow the format from the current matrix. */
    Matrix matrix(this->dimensions.first, this->dimensions.second);
    for (unsigned int mi = 0; mi < matrix.getNRows(); mi++)
    {
        for (unsigned int ni = 0; ni < matrix.getNColumns(); ni++)
        {
            matrix.rawMatrix[mi][ni] = this->rawMatrix[mi][ni] - M.rawMatrix[mi][ni];
        }
    }

    return matrix;
}

Matrix Matrix::operator * (Matrix &M)
{
    /* Follow the format from the current matrix. */
    Matrix matrix(this->dimensions.first, M.dimensions.second, 0);
    for (unsigned int i = 0; i < matrix.getNRows(); i++)
    {
        for (unsigned int j = 0; j < matrix.getNColumns(); j++)
        {
            for (unsigned int k = 0; k < this->getNColumns(); k++)
            {
                matrix.rawMatrix[i][j] += this->rawMatrix[i][k] * M.rawMatrix[k][j];
            }
        }
    }

    return matrix;
}
