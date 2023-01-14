#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <functional>

void heatEvolveExplicitEulerPeriodic(std::function<double(double)> u_initial, double a = -1, double b = 1, unsigned int J = 20, unsigned int T = 10000, unsigned int M = 10000000);

#endif  // solvers.hpp
