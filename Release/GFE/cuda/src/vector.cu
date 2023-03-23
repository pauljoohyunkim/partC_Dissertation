#include "solver.cuh"
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

__device__ double operator * (Vector v1, Vector v2)
{
    double xp = v1.x * v2.x;
    double yp = v1.y * v2.y;
    double zp = v1.z * v2.z;
    return xp + yp + zp;
}

__device__ Vector operator ^ (Vector v1, Vector v2)
{
    double x = v1.y * v2.z - v1.z * v2.y;
    double y = v1.z * v2.x - v1.x * v2.z;
    double z = v1.x * v2.y - v1.y * v2.x;

    Vector v { x, y, z };

    return v;
}

__device__ double norm(Vector v)
{
    return norm3d(v.x, v.y, v.z);
}

__device__ Vector vectorFromTensor(double* dev_blocks, int i, unsigned int N)
{
    Vector v {
        componentAccess(dev_blocks, i, 0, N),
        componentAccess(dev_blocks, i, 1, N),
        componentAccess(dev_blocks, i, 2, N)
    };
    return v;
}
