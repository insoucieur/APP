#include "Matrix2x2.hpp"
#include <cmath>
#include <iostream>

// Override default constructor
// initialises all entries of matrix to 0
Matrix2x2::Matrix2x2()
{
    double A[2][2] = {{val00, val01}, {val10, val11}};
    val00 = 0.0;
    val01 = 0.0;
    val10 = 0.0;
    val11 = 0.0;

}

//an overridden copy constructor
Matrix2x2::Matrix2x2(const Matrix2x2& other)
{
    val00 = other.val00;
    val01 = other.val01;
    val10 = other.val10;
    val11 = other.val11;
}
//a constructor that specifies the four entries of the matrix
Matrix2x2::Matrix2x2(double a, double b, double c, double d)
{
    val00 = a;
    val01 = b;
    val10 = c;
    val11 = d;
}
//method (function) that returns determinant of matrix
double Matrix2x2::CalcDeterminant() const
{
    double det = (val00*val11) - (val01*val10);
    return det;
}
//method that returns inverse of matrix it exists i.e. det(A) != 0
Matrix2x2 Matrix2x2::CalcInverse() const
{
    double determinant = CalcDeterminant();
    if (determinant != 0)
    {
        double invDet = (double(1.0)/determinant);
        double Ia = invDet*val11;
        double Ib = -invDet*val01;
        double Ic = -invDet*val10;
        double Id = invDet*val00;
        Matrix2x2 inv(Ia, Ib, Ic, Id);
        return inv;
    }
}
//overload assignment operator and allows A = B for class A and B
Matrix2x2& Matrix2x2::operator=(const Matrix2x2& z)
{
    val00 = z.val00;
    val01 = z.val01;
    val10 = z.val10;
    val11 = z.val11;
    return *this;
}
//overload unary subtraction operator so can write A = - B for class A and B
Matrix2x2 Matrix2x2::operator-() const
{
    Matrix2x2 B;
    B.val00 = -val00;
    B.val01 = -val01;
    B.val10 = -val10;
    B.val11 = -val11;
    return B;
}
//overload binary addition and subtraction operators so can write A = B + C
//or A = B -C for class A, B, C
Matrix2x2 Matrix2x2::operator+(const Matrix2x2& z) const
{
    Matrix2x2 C;
    C.val00 = val00 + z.val00;
    C.val01 = val01 + z.val01;
    C.val10 = val10 + z.val10;
    C.val11 = val11 + z.val11;
    return C;

}
Matrix2x2 Matrix2x2::operator-(const Matrix2x2& z) const
{
    Matrix2x2 C;
    C.val00 = val00 - z.val00;
    C.val01 = val01 - z.val01;
    C.val10 = val10 - z.val10;
    C.val11 = val11 - z.val11;
    return C;
}

//method that multiples a matrix by a specific double float variable
void Matrix2x2::MultScalar(double scalar)
{
    val00 = scalar*val00;
    val01 = scalar*val01;
    val10 = scalar*val10;
    val11 = scalar*val11;
}
