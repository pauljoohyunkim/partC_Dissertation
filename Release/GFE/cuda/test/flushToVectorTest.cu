#include <vector>
#include "../src/solver.cuh"

int main()
{
    std::vector<double> x { 1, 2, 3 };
    std::vector<double> y { 4, 5, 6 };
    std::vector<double> z { 7, 8, 9 };

    CurveTensor T { x, y, z };

    /* Writing to new vector */
    std::vector<double> ax;
    std::vector<double> ay;
    std::vector<double> az;
    ax.resize(3);
    ay.resize(3);
    az.resize(3);
    
    tensorBlockFlush(T, ax, ay, az);

    return 0;
}
