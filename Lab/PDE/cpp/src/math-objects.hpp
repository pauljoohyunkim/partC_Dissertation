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

    private:
        /* Matrix Characteristics */
        std::pair<unsigned int, unsigned int> dimensions {1, 1};
        std::vector<std::vector<double>> rawMatrix;


};


#endif  // MATH_OBJECTS_HPP
