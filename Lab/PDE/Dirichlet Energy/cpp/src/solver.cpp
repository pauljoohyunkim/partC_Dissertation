#include "solver.hpp"

Solver::Solver(std::function<double(double)> au_initial, double aa, double ab, unsigned int aJ, unsigned int aT, unsigned int aM)
{
    u_initial = au_initial;
    a = aa;
    b = ab;
    J = aJ;
    T = aT;
    M = aM;
}
