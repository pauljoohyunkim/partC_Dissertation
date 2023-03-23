#include <vector>
#include <cstdio>
#include <cmath>
#include "../src/tpe.cuh"
#include "../src/solver.cuh"

__global__ void kernel(double* dev_blocks, unsigned int N)
{
   double xi;
   double eta;
   Vector dxi;
   Vector deta;
   int p = 0;
   int q = 2;
   int r = 0;

   kjk(dev_blocks, p, q, r, N, xi, eta, dxi, deta);

   printf("xi: %f\n", xi);
   printf("eta: %f\n", eta);
   printf("dxi: %f, %f, %f\n", dxi.x, dxi.y, dxi.z);
   printf("deta: %f, %f, %f\n", deta.x, deta.y, deta.z);
}

int main()
{
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};
    const int res = 40;

    for (auto i = 0; i < res; i++)
    {
        double theta = 2 * M_PI / res * i;
        x.push_back(cos(theta));
        y.push_back(sin(2*theta));
        z.push_back(0.2 * sin(theta));
    }
    CurveTensor T { x, y, z };

    kernel<<<1,1>>>(T.dev_blocks, T.N);

    return 0;
}
