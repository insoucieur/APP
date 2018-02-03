#include "6.1.2.hpp"
#include <cmath>
#include <iostream>

int main()
{
    ComplexNumber z0(0.0, 0.0);
    ComplexNumber** A = new ComplexNumber* [3];
    ComplexNumber** res = new ComplexNumber* [3];
    for (int i = 0; i< 3; i++ )
    {
        A[i] = new ComplexNumber [3];
        res[i] = new ComplexNumber [3];
    }
    ComplexNumber z1(2.0,3.0);
    ComplexNumber z2(4.0,5.0);
    ComplexNumber z3(6.0,7.0);
    A[0][0] = z1;
    A[1][0] = z0;
    A[2][0] = z0;
    A[0][2] = z0;
    A[1][2] = z0;
    A[2][2] = z3;
    A[0][1] = z0;
    A[1][1] = z2;
    A[2][1] = z0;
    std::cout << "A00 = " << A[0][0] << "\n";
    std::cout << "A11 = " << A[1][1] << "\n";
    std::cout << "A22 = " << A[2][2] << "\n";
    CalculateExponential(A, 20, res);
    std::cout << "res00 = " << res[0][0] << "\n";
    std::cout << "res11 = " << res[1][1] << "\n";
    std::cout << "res10 = " << res[1][0] << "\n";
    for (int i = 0; i<3; i++)
        {
            delete[] A[i]; //delete array within matrix
            delete[] res[i];
        }
        delete[] A;
        delete[] res;
}
