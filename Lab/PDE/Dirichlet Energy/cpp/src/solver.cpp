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

void Solver::setScheme(std::function<void(std::function<double(double)>, double, double, unsigned int, unsigned int, unsigned int, double**&)> aScheme)
{
    schemeFun = aScheme;
}

void Solver::solve()
{
    if(!qSolved)
    {
        qSolved = true;
        Solver::allocate();
        schemeFun(u_initial, a, b, J, T, M, u);
    }
}

void Solver::plotSolution()
{
    if (!qSolved)
    {
        Solver::solve();
    }
}

void Solver::allocate()
{
    if(!qAllocated)
    {
        u = new double* [M + 1];
        for (auto m = 0; m <= M; m++)
        {
            u[m] = new double [J + 1];
        }
        qAllocated = false;
    }
}

void Solver::free()
{
    if(!qAllocated)
    {
        for (auto m = 0; m <= M; m++)
        {
            delete [] (u[m]);
        }
        delete [] u;
    }
}
