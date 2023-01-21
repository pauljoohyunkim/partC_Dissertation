#include "math-objects.hpp"
#include <stdexcept>

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
    if (this->dimensions != M.dimensions)
    {
        throw std::length_error("Matrix addition for given matrices is not defined.");
    }

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
    if (this->dimensions != M.dimensions)
    {
        throw std::length_error("Matrix subtraction for given matrices is not defined.");
    }

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
    if (this->dimensions.second != M.getNRows())
    {
        throw std::length_error("Matrix multiplication for given matrices is not defined.");
    }

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

std::vector<double>& Matrix::operator [](unsigned int i)
{
    return this->rawMatrix[i];
}

///* Vector Class Constructor */
Vector::Vector(unsigned int n, double initval) : Matrix::Matrix(n, 1, initval)
{
}

Vector::Vector(std::vector<double> &stdvector) : Matrix::Matrix(stdvector.size(), 1)
{
    for (unsigned int i = 0; i < stdvector.size(); i++)
    {
        rawMatrix[i][0] = stdvector[i];
    }
}

double& Vector::operator [] (unsigned int i)
{
    return this->rawMatrix[i][0];
}

double Vector::operator % (Vector& v)
{
    double x {0};
    if (this->dimensions.first != v.dimensions.first)
    {
        throw std::length_error("Scalar product not defined on vectors of different lengths.");
    }

    /* Sum */
    for (unsigned int i = 0; i < v.dimensions.first; i++)
    {
        x += this->rawMatrix[i][0] * v[i];
    }

    return x;
}

Vector3D::Vector3D(double initval) : Vector::Vector(3, initval)
{
}

Vector3D::Vector3D(double v1, double v2, double v3) : Vector::Vector(3)
{
    rawMatrix[0][0] = v1;
    rawMatrix[1][0] = v2;
    rawMatrix[2][0] = v3;
}

Vector3D::Vector3D(std::vector<double> &stdvector) : Vector::Vector(stdvector)
{
}

Vector3D Vector3D::operator ^ (Vector3D &v)
{
    Vector3D cross(0);
    cross[0] = this->rawMatrix[1][0] * v[2] - this->rawMatrix[2][0] * v[1];
    cross[1] = this->rawMatrix[2][0] * v[0] - this->rawMatrix[0][0] * v[2];
    cross[2] = this->rawMatrix[0][0] * v[1] - this->rawMatrix[1][0] * v[0];

    return cross;

}

/* Functions */
/* Vectorization
 * vectorize function takes a Matrix (VARIABLE) and flattens it to a vertical vector.
 * vectorize 3d function takes a Matrix (VARIABLE) and takes the (0,0), (1,0), (2,0) entries
 * to create a 3D vertical vector.
 * */
Vector vectorize(Matrix &M)
{
    auto nRow = M.getNRows();
    auto nColumn = M.getNColumns();

    
    Vector v(M.getNRows() * M.getNColumns());
    for (unsigned int mi = 0; mi < nRow; mi++)
    {
        for (unsigned int ni = 0; ni < nColumn; ni++)
        {
            v[mi + nColumn * ni] = M[mi][ni];
        }
    }

    return v;
}


Vector3D vectorize3d(Matrix &v)
{
    Vector3D v3(v[0][0], v[1][0], v[2][0]);

    return v3;
}


