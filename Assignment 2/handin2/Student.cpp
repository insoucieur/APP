#include "Student.hpp"
#include <cmath>
#include <iostream>
#include <string>

Student::Student()
{
    name = "unspecified";
    library_fines = 0.0;
    tuition_fees = 0.0;
}
Student::Student(std::string name, double fines, double fees)
{
    this->name = name;
    tuition_fees = fees;
    library_fines = fines;
}

double Student::MoneyOwed() const
{
    return (library_fines + tuition_fees);
}

//method allows user to set library fines only to non -ve values
void Student::SetLibraryFines(double amount)
{
    library_fines = fabs(amount);
}
//another method to access the private variable libraryfines
double Student::GetLibraryFines() const
{
    return library_fines;
}


