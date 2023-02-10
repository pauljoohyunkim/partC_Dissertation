#ifndef __CURVE_HPP__
#define __CURVE_HPP__

#include <vector>

class FourierCurve
{
    public:
        FourierCurve(std::vector<double>& axa, std::vector<double>& axb, std::vector<double>& aya, std::vector<double>& ayb, std::vector<double>& aza, std::vector<double>& azb);
        ~FourierCurve();

        /* Load data onto GPU VRAM */
        void cudafy();
        
        unsigned int J { 0 };
        std::vector<double> xa;
        std::vector<double> xb;
        std::vector<double> ya;
        std::vector<double> yb;
        std::vector<double> za;
        std::vector<double> zb;

        /* Pointer to memory on GPU */
        double* dev_coefficients;
        double* dev_scratch_pad;

    private:
        bool dev_coefficient_allocated { false };
        bool dev_scratch_pad_allocated { false };
};

/* Primitive DEBUG functions */
__global__ void printCoefficientsPartiallyDEBUG(double* devcoeffs, int index);

#endif
