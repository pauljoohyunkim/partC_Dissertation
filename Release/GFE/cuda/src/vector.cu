#include "vector.cuh"


__device__ Vector::Vector()
{
    x = 0;
    y = 0;
    z = 0;
}

__device__ Vector::Vector(double ax, double ay, double az)
{
    x = ax;
    y = ay;
    z = az;
}

__device__ double Vector::norm()
{
    return norm3d(x, y, z);
}

__device__ Vector operator + (Vector v1, Vector v2)
{
    Vector v { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    return v;
}

__device__ Vector operator - (Vector v1, Vector v2)
{
    Vector v { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    return v;
}
