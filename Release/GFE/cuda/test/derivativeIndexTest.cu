#include "../src/solver.cuh"
#include "../src/tpe.cuh"
#include <cstdio>

__global__ void kernel(int* dev_derivative_indices, unsigned int N)
{
    fillDerivativeIndex(dev_derivative_indices, 1, N);

    for (int i = 0; i < 4 * ((int) N - 3); i++)
    {
        printf("(%d, %d)\n", dev_derivative_indices[2 * i], dev_derivative_indices[2 * i + 1]);
    }
}

int main()
{
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};
    const int res = 6;

    for (auto i = 0; i < res; i++)
    {
        double theta = 2 * M_PI / res * i;
        x.push_back(cos(theta));
        y.push_back(sin(2*theta));
        z.push_back(0.2 * sin(theta));
    }
    CurveTensor T { x, y, z };

    kernel<<<1,1>>>(T.dev_derivative_indices, T.N);

    cudaDeviceSynchronize();
    return 0;
}
