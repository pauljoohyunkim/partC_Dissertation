#ifndef MATH_OBJECTS_HPP
#define MATH_OBJECTS_HPP

#include <utility>
#include <vector>

class Matrix
{
    public:
        /* Matrix Constructor 
         * First constructor generates a matrix with all entries equal to initval
         * Second constructor generates a matrix from double pointer pMatrix
         * */

        Matrix(unsigned int m, unsigned int n, double initval = 0);
        Matrix(unsigned int m, unsigned int n, double** &pMatrix);

        /* Class Functions */
        unsigned int getNRows();
        unsigned int getNColumns();
        std::pair<unsigned int, unsigned int> getNRowsColumns();

        /* Operators */
        Matrix operator + (Matrix &M);
        Matrix operator - (Matrix &M);
        Matrix operator * (Matrix &M);
        Matrix operator * (double lambda);
        std::vector<double>& operator [](unsigned int i);

    protected:
        /* Matrix Characteristics */
        std::pair<unsigned int, unsigned int> dimensions {1, 1};
        std::vector<std::vector<double>> rawMatrix;

};

class Vector: public Matrix
{
    public:
        /* Vector Constructor */
        Vector(unsigned int n, double initval = 0);
        Vector(std::vector<double> &stdvector);
        /* Operator */
        Vector operator + (Vector &v);
        Vector operator - (Vector &v);
        Vector operator * (double lambda);
        /* Fetch values instead of std::vector */
        double& operator [] (unsigned int i);
        /* Scalar Product */
        double operator % (Vector& v);
};

class Vector3D: public Vector
{
    public:
        /* Vector3D Constructor */
        Vector3D(double initval = 0);
        Vector3D(double v1, double v2, double v3);
        Vector3D(std::vector<double> &stdvector);

        /* Operator */
        Vector3D operator + (Vector3D &v);
        Vector3D operator - (Vector3D &v);
        Vector3D operator * (double lambda);
        Vector3D operator ^ (Vector3D& v);
};

/* Get column vector */
Vector columnvector(Matrix &M, unsigned int i);
Vector3D columnvector3D(Matrix &M, unsigned int i);

/* Matrix-Vector Multiplication */
Vector matvecmul(Matrix &M, Vector &v);
Vector3D matvecmul3D(Matrix &M, Vector3D &v);

/* Vectorization
 * vectorize function takes a Matrix (VARIABLE) and flattens it to a vertical vector.
 * vectorize 3d function takes a Matrix (VARIABLE) and takes the (0,0), (1,0), (2,0) entries
 * to create a 3D vertical vector.
 * */
Vector vectorize(Matrix &M);
Vector3D vectorize3d(Matrix &v);

/* Norm */
double l2norm(Vector& v);



#endif  // MATH_OBJECTS_HPP
