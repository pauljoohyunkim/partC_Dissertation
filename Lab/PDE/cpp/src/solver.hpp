#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <functional>
#include <iostream>
#include <vector>

#define DOUBLE_PRECISION 15

class Solver
{
    public:
        /* Constructor */
        Solver(std::function<double(double)> au_initial, double aa = -1, double ab = 1, unsigned int aJ = 20, double aT = 1, unsigned int aM = 1000);

        /* Destructor */
        ~Solver();

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

#endif  // solver.hpp
