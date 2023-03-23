#include "../src/tpe.cuh"
#include <cstdio>

__global__ void kernel()
{
    Vector T {-0.4933978725646011, 0.8547676525893596, 0.16103043015406};
    Vector p {0.8910065241883679, 0.8090169943749475, 0.09079809994790936};
    Vector q {0.9876983405951377, 0.3090169943749474, 0.03128689300804618};
    auto kalphabeta = kernelalphabeta(p, q, T, 3, 6);

    printf("%f\n", kalphabeta);
}

int main()
{
    kernel<<<1,1>>>();

    cudaDeviceSynchronize();

    return 0;
}
