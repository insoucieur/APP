#include "GraduateStudent.hpp"
#include <cmath>
#include <iostream>
#include <string>

//use polymorphism to write a method that calculate total money
//owed by a graduate student, who do not pay tuition fees
//will require method for calculating total money owed to be a virtual function
//of parent class
GraduateStudent::GraduateStudent()
{

}
GraduateStudent::GraduateStudent(std::string name, double fines, double fees, bool fullTime):Student(name, fines, fees)
{
    fullTime = fullTime;
}
double GraduateStudent::MoneyOwed() const
{
    return Student::GetLibraryFines(); //way to access private members in derived class?
}

