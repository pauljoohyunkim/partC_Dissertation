#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <functional>

class Solver
{
    public:
        /* Constructor */
        Solver(std::function<double(double)> au_initial, double aa = -1, double ab = 1, unsigned int aJ = 20, unsigned int aT = 1, unsigned int aM = 1000);

        /* Set scheme */
        void setScheme(std::function<void(std::function<double(double)>, double, double, unsigned int, unsigned int, unsigned int)> aScheme);

        /* Solve */
        void solve();

        /* Plot */
        void plotSolution();
        
        
    private:
        /* Interval of interest: (a, b) */
        double a;
        double b;
        unsigned int J;         /* The last spacial step */
        unsigned int T;         /* Time interval (0,T) */
        unsigned int M;         /* The last time step */
        double** u;             /* Entire Solution */
        bool qSolved {false};   /* Whether or not if solve was used before */
        std::function<double(double)> u_initial;    /* Initial Datum */
        /* Scheme function that takes u_initial, a, b, J, T, and M as parameters */
        std::function<void(std::function<double(double)>, double, double, unsigned int, unsigned int, unsigned int)> schemeFun;

};

#endif  // solver.hpp
