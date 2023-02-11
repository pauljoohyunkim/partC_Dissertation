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

        /* Pull data from GPU VRAM to host */
        void cudaFlush();

        
        unsigned int J { 0 };
        std::vector<double> xa;
        std::vector<double> xb;
        std::vector<double> ya;
        std::vector<double> yb;
        std::vector<double> za;
        std::vector<double> zb;
        std::vector<double> coeff_differential;
        unsigned int resolution { DEFAULT_RESOLUTION };
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

        /* Pointer to memory on GPU */
        /* dev_coefficients is a concatenation of xa, xb, ya, yb, za, zb */
        /* dev_trig_val_table is a table of precomputed cosine and sine values. */
        double* dev_coefficients;
        double* dev_differential_coefficients;
        double* dev_cos_table;
        double* dev_sin_table;
        double* dev_x;
        double* dev_y;
        double* dev_z;

    private:
        bool dev_coefficient_allocated { false };
        bool dev_differential_coefficient_allocated { false };
        bool dev_trig_val_table_allocated { false };
        bool dev_curve_points_allocated { false };
};

/* Cross Product */
__device__ void cross(double x1, double x2, double x3, double y1, double y2, double y3, double& z1, double& z2, double& z3);
/* trig_table query
 * dev_table is either dev_cos_table or dev_sin_table */
__device__ double dev_trig_table_query(double* dev_table, unsigned int i, unsigned int k);
/* Fill curve position vectors from coefficients */
__device__ void fill_pos(double* dev_x, double* dev_y, double* dev_z, double* dev_coefficients, double* dev_cos_table, double* dev_sin_table, unsigned int resolution, unsigned int J);

/* Primitive DEBUG functions */
__global__ void printCoefficientsPartiallyDEBUG(double* device_float_value);
__global__ void crossDEBUG(double x1, double x2, double x3, double y1, double y2, double y3);
__global__ void queryDEBUG(double* dev_table, int i, int k, unsigned int J);
__global__ void fillDEBUG(double* dev_x, double* dev_y, double* dev_z, double* dev_coefficients, double* dev_cos_table, double* dev_sin_table, unsigned int resolution, unsigned int J);

#endif
