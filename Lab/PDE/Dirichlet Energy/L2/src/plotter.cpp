#include "plotter.hpp"
#include <matplot/matplot.h>

void plot(std::vector<double> x, std::vector<double> y)
{
    matplot::plot(x,y);
}
