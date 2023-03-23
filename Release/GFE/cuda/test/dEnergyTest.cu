#include "../src/solver.cuh"
#include "../src/tpe.cuh"
#include "../src/vector.cuh"
#include <vector>

int main()
{
    const unsigned int resolution = 30;
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> z {};
    for (auto i = 0; i < resolution; i++)
    {
        double theta = 2 * M_PI / resolution * i;
        x.push_back(cos(theta));
        y.push_back(sin(2*theta));
        z.push_back(0.2 * sin(theta));
    }

    /* Curve and Gradient */
    CurveTensor T { x, y, z };
    CurveTensor dT { resolution };

    /* Derivative Index */
    ScratchPad<int> s { resolution, 8 * (resolution - 3) };

    cuDEnergy<<<resolution, 1>>>(T.dev_blocks, dT.dev_blocks, s.scratchpads, resolution, 3, 6);

    /* Get Differential into RAM */
    double derivatives_in_memory[3 * resolution];
    tensorBlockFlush(dT, derivatives_in_memory, resolution);


    return 0;

    
}
