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
        double& operator [](unsigned int i);
};

class Vector3D: public Vector
{
    public:
        /* Vector3D Constructor */
        Vector3D(double initval = 0);
        Vector3D(double v1, double v2, double v3);
        Vector3D(std::vector<double> &stdvector);

        /* Operator */
        //Vector3D operator ^ (Vector3D& v);
};


#endif  // MATH_OBJECTS_HPP
