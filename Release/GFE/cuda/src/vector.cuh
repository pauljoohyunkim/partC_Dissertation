#ifndef __VECTOR_CUH__
#define __VECTOR_CUH__

struct Vector
{
    __device__ Vector();
    __device__ Vector(double ax, double ay, double az);

    /* Norm */
    __device__ double norm();
    
    double x;
    double y;
    double z;
};

__device__ Vector operator + (Vector v1, Vector v2);
__device__ Vector operator - (Vector v1, Vector v2);
/* Dot Product */
__device__ double operator * (Vector v1, Vector v2);
__device__ Vector operator ^ (Vector v1, Vector v2);

#endif  // __VECTOR_CUH__
