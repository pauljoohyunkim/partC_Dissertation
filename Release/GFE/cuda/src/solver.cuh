#ifndef __SOLVER_CUH__
#define __SOLVER_CUH__

class CurveTensor
{
    public:
        CurveTensor(unsigned int aN);
        ~CurveTensor();

        unsigned int N { 0 };

        double* dev_blocks;
};

void tensorBlockLoad(CurveTensor& Gammabf, double* blocks, unsigned int N);
void tensorBlockFlush(CurveTensor& Gammabf, double* blocks, unsigned int N);
void tensorAdd(CurveTensor& t1, CurveTensor& t2);

#endif
