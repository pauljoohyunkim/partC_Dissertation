#include "../src/math-objects.hpp"

int main()
{
     Matrix A(3,2);

     double** b = new double* [3];
     b[0] = new double [2];
     b[1] = new double [2];
     b[2] = new double [2];
     

     b[0][0] = 3;
     b[0][1] = 2;
     b[1][0] = 6;
     b[1][1] = 4.2;
     b[2][0] = 1;

     Matrix B(3, 2, b);

     return 0;
}
