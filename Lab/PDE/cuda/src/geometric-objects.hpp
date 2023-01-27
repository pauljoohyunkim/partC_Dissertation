#ifndef GEOMETRIC_OBJECTS_H
#define GEOMETRIC_OBJECTS_H

#include <vector>

class cuCurve
{
    public:
        /* Constructor */
        cuCurve(unsigned int aJ);
        cuCurve(std::vector<double> &aX, std::vector<double> &aY, std::vector<double> &aZ);

        /* Deconstructor */
        ~cuCurve();

        /* Variables */
        unsigned int J;

        /* In order to use device functions, one needs to turn the array into CUDA array by
         * invoking cudafy in the beginning */
        void cudafy();
        /* Flush progress to x, y, z from GPU */
        void flushFromDevice();
        /* CUDA Functions */
        __device__ double& getX(int i);
        __device__ double& getY(int i);
        __device__ double& getZ(int i);
        double* dev_x {};
        double* dev_y {};
        double* dev_z {};

        /* Non-cuda functions, but requires cudafy to have been called. */
        double getValFromDevice(int coordnum, int i);
    protected:
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        bool dev_x_allocated { false };
        bool dev_y_allocated { false };
        bool dev_z_allocated { false };
};


#endif  // geometric-objects.h
