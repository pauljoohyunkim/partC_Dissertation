#include "solver.hpp"
#include "params.hpp"
#include <iostream>
#include <fstream>
#include <matplot/matplot.h>
#include <iomanip>

Solver::Solver(std::function<double(double)> au_initial, double aa, double ab, unsigned int aJ, unsigned int aT, unsigned int aM)
{
    u_initial = au_initial;
    a = aa;
    b = ab;
    J = aJ;
    T = aT;
    M = aM;
}

Solver::~Solver()
{
    free();
    std::cout << "Solver Destructed" << std::endl;
}

void Solver::setScheme(std::function<void(std::function<double(double)>, double, double, unsigned int, unsigned int, unsigned int, double**&, std::vector<double>&)> aScheme)
{
    schemeFun = aScheme;
}

void Solver::solve()
{
    if(!qSolved)
    {
        qSolved = true;
        Solver::allocate();
        schemeFun(u_initial, a, b, J, T, M, u, x);
    }
}

/* DEPRECATED!!!!!!
 * Use exportSolution, then the python program json_plot.py to get a visual solution! */
void Solver::plotSolution()
{
    std::vector<double> y {};
    double deltaX = (b - a) / J;
    if (!qAllocated || !qSolved)
    {
        Solver::solve();
    }

    /* Preallocation of y */
    for (unsigned int j = 0; j <= J; j++)
    {
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

    /* Exporting data in json file 
     * First entry is the list of x values.
     * Second entry is the 2D array of solution in the form u[m][j]
     * */

    exportFile << "[";

    exportFile << "[";
    for(unsigned int j = 0; j <= J; j++)
    {
        exportFile << std::setprecision(DOUBLE_PRECISION) << x[j];
        if(j != J)
        {
            exportFile << ",";
        }
    }
    exportFile << "]";
    


    exportFile << "," << std::endl;

    exportFile << "[";
    for(unsigned int m = 0; m <= M; m++)
    {
        exportFile << "[";
        for(unsigned int j = 0; j <= J; j++)
        {
            if(j == J)
            {
                exportFile << std::setprecision(DOUBLE_PRECISION) << u[m][j];
            }
            else
            {
                exportFile << std::setprecision(DOUBLE_PRECISION) <<  u[m][j] << ",";
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
