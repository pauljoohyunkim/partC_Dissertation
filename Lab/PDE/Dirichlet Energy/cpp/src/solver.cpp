#include "solver.hpp"
#include "params.hpp"
#include <iostream>
#include <fstream>
#include <matplot/matplot.h>

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
    std::vector<double> x {};
    std::vector<double> y {};
    double deltaX = (b - a) / J;
    if (!qAllocated || !qSolved)
    {
        Solver::solve();
    }

    /* x values and preallocation of y */
    for (unsigned int j = 0; j <= J; j++)
    {
        x.push_back(a + (double) j * deltaX);
        y.push_back(0);
    }

    matplot::hold(matplot::on);
    /* Load each array to vector, then plot */
    for (unsigned int m = 0; m <= M; m++)
    {
        y.assign(u[m], u[m] + J + 1);

        matplot::plot(x, y);
        matplot::xlim({XLIM_L, XLIM_R});
        matplot::ylim({YLIM_D, YLIM_U});
    }

    /* Deallocate */
    Solver::free();
}

void Solver::exportSolution(std::string filename)
{
    std::ofstream exportFile(filename);

    /* Exporting data in json file */
    exportFile << "[";
    for(unsigned int m = 0; m <= M; m++)
    {
        exportFile << "[";
        for(unsigned int j = 0; j <= J; j++)
        {
            if(j == J)
            {
                exportFile << u[m][j];
            }
            else
            {
                exportFile << u[m][j] << ",";
            }
        }
        if (m == M)
        {
            exportFile << "]";
        }
        else
        {
            exportFile << "]," << std::endl;
        }
    }
    exportFile << "]";

    exportFile.close();

    /* Deallocate */
    Solver::free();
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
        qAllocated = true;
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
    qAllocated = false;
}
