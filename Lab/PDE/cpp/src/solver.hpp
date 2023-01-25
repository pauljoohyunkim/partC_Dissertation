#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <functional>
#include <iostream>
#include <vector>
#include "math-objects.hpp"
#include "geometric-objects.hpp"

#define DOUBLE_PRECISION 15

/* PDE from Dirichlet energy Solver */
class Solver1D
{
    public:
        /* Constructor */
        Solver1D(std::function<double(double)> au_initial, double aa = -1, double ab = 1, unsigned int aJ = 20, double aT = 1, unsigned int aM = 1000);

        /* Destructor */
        ~Solver1D();

        /* Set scheme */
        void setScheme(std::function<void(std::function<double(double)>, double, double, unsigned int, double, unsigned int, double**&, std::vector<double>&)> aScheme);

        /* Solve */
        void solve();

        /* Export */
        void exportSolution(std::string filename, unsigned int timeskip = 1);
        
        
        /* Mesh Points */
        std::vector<double> x;
    private:
        /* Interval of interest: (a, b) */
        double a;
        double b;
        unsigned int J;         /* The last spacial step */
        double T;               /* Time interval (0,T) */
        unsigned int M;         /* The last time step */
        double** u;             /* Entire Solution */
        bool qAllocated {false};    /* Whether or not if u is allocated or not */
        bool qSolved {false};   /* Whether or not if solve was used before */
        std::function<double(double)> u_initial;    /* Initial Datum */
        /* Scheme function that takes u_initial, a, b, J, T, and M as parameters */
        std::function<void(std::function<double(double)>, double, double, unsigned int, double, unsigned int, double**&, std::vector<double>&)> schemeFun;

        /* Allocator for u */
        void allocate();

        /* Free for u */
        void free();
};

class SolverCurveRepulsion
{
    public:
        /* SolverCurveRepulsion class Constructor */
        /* For kernel function k_{ij}, arguments are x_i, x_{i+1}, x_j, x_{j+1}, T_{i}, and the Discrete Kernel object itself */
        SolverCurveRepulsion(double aAlpha, double aBeta, std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, SolverCurveRepulsion&)> aKernelFunction);

        /* Discretized Energy */
        double energy(Curve C);
        double energyDifferential(Curve C, Vector3D v, int index);



        /* k_{i,j} */
        std::function<double(Vector3D&, Vector3D&, Vector3D&, Vector3D&, Vector3D&, SolverCurveRepulsion&)> kernelFunction;
        /* k_\beta^\alpha (Actual) */
        std::function<double(Vector3D&, Vector3D&, Vector3D&)> kernelalphabeta;

    private:
        double alpha { 2 };
        double beta { 4 };
};

#endif  // solver.hpp
