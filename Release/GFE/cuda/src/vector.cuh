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
__device__ Vector operator - (Vector v1);
/* Scalar Multiplication */
__device__ Vector operator * (Vector v1, double lambda);
__device__ Vector operator * (double lambda, Vector v1);
/* Dot Product */
__device__ double operator * (Vector v1, Vector v2);
/* Cross Product */
__device__ Vector operator ^ (Vector v1, Vector v2);

__device__ double norm(Vector v);

__device__ Vector vectorFromTensor(double* dev_blocks, int i, unsigned int N);

#endif  // __VECTOR_CUH__
