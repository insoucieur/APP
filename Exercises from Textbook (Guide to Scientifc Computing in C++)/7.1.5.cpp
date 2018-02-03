#include "PhdStudent.hpp"
#include <cmath>
#include <iostream>
#include <string>
PhdStudent::PhdStudent(std::string name, double fines, double fees, bool fullTime):GraduateStudent(name,fines,fees,fullTime)
{

}
double PhdStudent::MoneyOwed() const
{
    return 0.0;
}
