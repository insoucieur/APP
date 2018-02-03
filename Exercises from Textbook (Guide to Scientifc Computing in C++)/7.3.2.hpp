#include "7.3.hpp"

#ifndef FORWARDEULERSOLVERDEF
#define FORWARDEULERSOLVERDEF

class ForwardEulerSolver: public AbstractOdeSolver
{
public:
    double RightHandSide (double y, double t);
    double SolveEquation();
};

#endif
