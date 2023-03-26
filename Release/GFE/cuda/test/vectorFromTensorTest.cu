#include "../src/vector.cuh"
#include "../src/solver.cuh"
#include <cstdio>
#include <vector>


__global__ void kernel(double* blocks, unsigned int N)
{
    auto v = vectorFromTensor(blocks, 2, N);

    printf("%f, %f, %f\n", v.x, v.y, v.z);
}

int main()
{
    std::vector x { 1.0, 4.2 };
    std::vector y { 2.0, 5.2 };
    std::vector z { 3.0, 2.2 };
    CurveTensor T { x, y, z };

    kernel<<<1,1>>>(T.dev_blocks, T.N);


    cudaDeviceSynchronize();
    return 0;
}
