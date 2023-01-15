#ifndef SOLVER_HPP
#define SOLVER_HPP

class Solver
{
    public:
        Solver(std::function<double(double)> au_initial, double aa = -1, double ab = 1, unsigned int aJ = 20, unsigned int aT = 1, unsigned int aM = 1000);
        
    private:
        /* Interval of interest: (a, b) */
        double a;
        double b;
        unsigned int J;         /* The last spacial step */
        unsigned int T;         /* Time interval (0,T) */
        unsigned int M;         /* The last time step */
        std::function<double(double)> u_initial;    /* Initial Datum */

};

#endif  // solver.hpp
