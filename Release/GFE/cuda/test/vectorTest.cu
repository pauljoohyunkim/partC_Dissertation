#include <cstdio>
#include "../src/vector.cuh"

__global__ void kernel()
{
    Vector v1 { 0, 1, 2 };
    Vector v2 { 1, -2, 5 };
    auto vp = v1 + v2;
    auto vn = v1 - v2;

    printf("v1: %f, %f, %f\n", v1.x, v1.y, v1.z);
    printf("v2: %f, %f, %f\n", v2.x, v2.y, v2.z);
    printf("sum: %f, %f, %f\n", vp.x, vp.y, vp.z);
    printf("diff: %f, %f, %f\n", vn.x, vn.y, vn.z);
}

int main()
{
    kernel<<<1,1>>>();
    
    cudaDeviceSynchronize();
    return 0;
}
