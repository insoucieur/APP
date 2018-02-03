#include <"7.3.hpp">

#ifndef RUNGEKUTTFOURTHORDERADEF
#define RUNGEKUTTFOURTHORDERADEF

class RungeKuttaFourthOrder: public AbstractOdeSolver
{
public:
    double RightHandSide (double y, double t);
    double SolveEquation();
};

#endif
