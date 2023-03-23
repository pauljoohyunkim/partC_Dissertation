#include <vector>
#include <cstdio>
#include <cmath>
#include <iostream>
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
   int k = 0;

   kjk(dev_blocks, p, q, r, N, xi, eta, dxi, deta);

   printf("kjk\n");
   printf("xi: %f\n", xi);
   printf("eta: %f\n", eta);
   printf("dxi: %f, %f, %f\n", dxi.x, dxi.y, dxi.z);
   printf("deta: %f, %f, %f\n", deta.x, deta.y, deta.z);

   p = 1;
   q = 4;
   r = 0;
   k = 0;

   ijk(dev_blocks, p, q, r, N, xi, eta, dxi, deta);
   
   printf("ijk\n");
   printf("xi: %f\n", xi);
   printf("eta: %f\n", eta);
   printf("dxi: %f, %f, %f\n", dxi.x, dxi.y, dxi.z);
   printf("deta: %f, %f, %f\n", deta.x, deta.y, deta.z);

   p = 0;
   q = 4;
   r = 0;
   k = 1;
    
   km1jkm1(dev_blocks, p, q, r, N, xi, eta, dxi, deta);
   printf("km1jkm1\n");
   printf("xi: %f\n", xi);
   printf("eta: %f\n", eta);
   printf("dxi: %f, %f, %f\n", dxi.x, dxi.y, dxi.z);
   printf("deta: %f, %f, %f\n", deta.x, deta.y, deta.z);

   p = 1;
   q = 4;
   r = 0;
   k = 1;

   kjkm1(dev_blocks, p, q, r, N, xi, eta, dxi, deta);
   printf("kjkm1\n");
   printf("xi: %f\n", xi);
   printf("eta: %f\n", eta);
   printf("dxi: %f, %f, %f\n", dxi.x, dxi.y, dxi.z);
   printf("deta: %f, %f, %f\n", deta.x, deta.y, deta.z);

   p = 3;
   q = 1;
   r = 4;
   k = 1;

   ikj(dev_blocks, p, q, r, N, xi, eta, dxi, deta);
   printf("ikj\n");
   printf("xi: %f\n", xi);
   printf("eta: %f\n", eta);
   printf("dxi: %f, %f, %f\n", dxi.x, dxi.y, dxi.z);
   printf("deta: %f, %f, %f\n", deta.x, deta.y, deta.z);

   p = 3;
   q = 1;
   r = 3;
   k = 1;

   Vector res { };
   dkalphabeta(dev_blocks, p, q, r, k, N, res, 3, 6);
   printf("%f, %f, %f\n", res.x, res.y, res.z);

}

__global__ void kernel2(double* dev_blocks, unsigned int N, int i, int j, int k)
{
    Vector res { };
    dkij(dev_blocks, i, j, k, N, res, 3, 6);
       printf("%f, %f, %f\n", res.x, res.y, res.z);
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
    cudaDeviceSynchronize();
    std::cout << "------------------------------" << std::endl;
    kernel2<<<1,1>>>(T.dev_blocks, T.N, 3, 1, 1);

    return 0;
}
