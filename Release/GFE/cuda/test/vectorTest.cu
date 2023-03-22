#include <cstdio>
#include "../src/vector.cuh"

__global__ void kernel()
{
    Vector v1 { 0, 1, 2 };
    Vector v2 { 1, -2, 5 };
    auto v = v1 + v2;

    printf("v1: %f, %f, %f\n", v1.x, v1.y, v1.z);
    printf("v2: %f, %f, %f\n", v2.x, v2.y, v2.z);
    printf("sum: %f, %f, %f\n", v.x, v.y, v.z);
}

int main()
{
    kernel<<<1,1>>>();
    
    cudaDeviceSynchronize();
    return 0;
}
