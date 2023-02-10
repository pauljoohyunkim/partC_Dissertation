#ifndef __CURVE_HPP__
#define __CURVE_HPP__

#include <vector>
#define DEFAULT_RESOLUTION 10

class FourierCurve
{
    public:
        FourierCurve(std::vector<double>& axa, std::vector<double>& axb, std::vector<double>& aya, std::vector<double>& ayb, std::vector<double>& aza, std::vector<double>& azb, unsigned int aresolution = DEFAULT_RESOLUTION);
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
        unsigned int resolution { DEFAULT_RESOLUTION };

        /* Pointer to memory on GPU */
        /* dev_coefficients is a concatenation of xa, xb, ya, yb, za, zb */
        /* dev_trig_val_table is a table of precomputed cosine and sine values. */
        double* dev_coefficients;
        double* dev_cos_table;
        double* dev_sin_table;

    private:
        bool dev_coefficient_allocated { false };
        bool dev_trig_val_table_allocated { false };
};

/* Primitive DEBUG functions */
__global__ void printCoefficientsPartiallyDEBUG(double* device_float_value);

#endif
