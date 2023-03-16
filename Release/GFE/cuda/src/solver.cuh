#ifndef __SOLVER_CUH__
#define __SOLVER_CUH__

#include <vector>

class CurveTensor
{
    public:
        CurveTensor(unsigned int aN);
        CurveTensor(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
        ~CurveTensor();

        unsigned int N { 0 };

        double* dev_blocks;
};

void tensorBlockLoad(CurveTensor& Gammabf, double* blocks, unsigned int N);
void tensorBlockFlush(CurveTensor& Gammabf, double* blocks, unsigned int N);
void tensorAdd(CurveTensor& t1, CurveTensor& t2);
void tensorSubtract(CurveTensor& t1, CurveTensor& t2);

#endif