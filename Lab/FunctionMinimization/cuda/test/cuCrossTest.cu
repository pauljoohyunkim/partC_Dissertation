#include "../src/curve.hpp"

int main()
{
    crossDEBUG<<<1,1>>>(2, -3, -1, 1, 4, -2);
    cudaDeviceSynchronize();

    return 0;
}
