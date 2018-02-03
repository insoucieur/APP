#include "CalculateExponential.hpp"
#include <cmath>
#include <iostream>
#include "ComplexNumber.hpp"
void CalculateExponential(ComplexNumber **A, int nMax, ComplexNumber **res){
ComplexNumber M[3][3];
     for (int i = 0; i< 3; i++ ){
        for (int j = 0; j< 3; j++ ){
            M[i][j] = A[i][j];
        }}
    ComplexNumber zZero;
    ComplexNumber zOne;
    ComplexNumber zTotal;
    ComplexNumber zPrev;
    ComplexNumber zCurrent;
    ComplexNumber zNext;
    ComplexNumber zSum;

for (int i = 0; i< 3; i++ ){
    for (int j = 0; j< 3; j++ ){
        //std::cout << "M00 = " << M[i][j] << "\n";
        zTotal = ComplexNumber (0.0, 0.0);
        int factorial = 1;
    if (i==j){
        zZero = ComplexNumber(1.0, 0.0);
        zOne = M[i][j];
        zPrev = M[i][j];
        zSum = ComplexNumber (0.0, 0.0);
    for(int n = 2; n <= nMax; n++)
{
                //std::cout << "Prev = " << zPrev << "\n";
                double a = RealPart(zPrev);
                double b = ImaginaryPart(zPrev);
                double c = RealPart(M[i][j])/double(n);
                double d = ImaginaryPart(M[i][j])/double(n);
                zNext = ComplexNumber(c,d);
                //std::cout << "Next = " << zNext << "\n";
                double Re = (a*c) - (b*d);
                double Im = (c*b) + (a*d);
                zCurrent = ComplexNumber(Re, Im);
                //std::cout << "Current = " << zCurrent << "\n";
                zSum = zSum + zCurrent;
                zPrev = zCurrent;
            }
            zTotal = zZero + zOne + zSum;
            //std::cout << "Total = " << zTotal << "\n";

        res[i][j] = zTotal;}
        else{
            res[i][j] = zTotal;
        }
        //std::cout << "Result = " << res[i][j] << "\n";
    }}
}

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
    std::cout << "res22 = " << res[2][2] << "\n";
    std::cout << "res10 = " << res[1][0] << "\n";
    for (int i = 0; i<3; i++)
        {
            delete[] A[i]; //delete array within matrix
            delete[] res[i];
        }
        delete[] A;
        delete[] res;
}



// Override default constructor
// Set real and imaginary parts to zero
ComplexNumber::ComplexNumber()
{
   mRealPart = 0.0;
   mImaginaryPart = 0.0;
}

// Constructor that sets complex number z=x+iy
ComplexNumber::ComplexNumber(double x, double y)
{
   mRealPart = x;
   mImaginaryPart = y;
}

// Method for computing the modulus of a
// complex number
double ComplexNumber::CalculateModulus() const
{
   return sqrt(mRealPart*mRealPart+
               mImaginaryPart*mImaginaryPart);
}

// Method for computing the argument of a
// complex number
double ComplexNumber::CalculateArgument() const
{
   return atan2(mImaginaryPart, mRealPart);
}

// Method for raising complex number to the power n
// using De Moivre's theorem - first complex
// number must be converted to polar form
ComplexNumber ComplexNumber::CalculatePower(double n) const
{
   double modulus = CalculateModulus();
   double argument = CalculateArgument();
   double mod_of_result = pow(modulus, n);
   double arg_of_result = argument*n;
   double real_part = mod_of_result*cos(arg_of_result);
   double imag_part = mod_of_result*sin(arg_of_result);
   ComplexNumber z(real_part, imag_part);
   return z;
}

// Overloading the = (assignment) operator
ComplexNumber& ComplexNumber::
               operator=(const ComplexNumber& z)
{
   mRealPart = z.mRealPart;
   mImaginaryPart = z.mImaginaryPart;
   return *this;
}

// Overloading the unary - operator
ComplexNumber ComplexNumber::operator-() const
{
   ComplexNumber w;
   w.mRealPart = -mRealPart;
   w.mImaginaryPart = -mImaginaryPart;
   return w;
}

// Overloading the binary + operator
ComplexNumber ComplexNumber::
              operator+(const ComplexNumber& z) const
{
   ComplexNumber w;
   w.mRealPart = mRealPart + z.mRealPart;
   w.mImaginaryPart = mImaginaryPart + z.mImaginaryPart;
   return w;
}

// Overloading the binary - operator
ComplexNumber ComplexNumber::
              operator-(const ComplexNumber& z) const
{
   ComplexNumber w;
   w.mRealPart = mRealPart - z.mRealPart;
   w.mImaginaryPart = mImaginaryPart - z.mImaginaryPart;
   return w;
}

// Overloading the insertion << operator
std::ostream& operator<<(std::ostream& output,
                         const ComplexNumber& z)
{
   // Format as "(a + bi)" or as "(a - bi)"
   output << "(" << z.mRealPart << " ";
   if (z.mImaginaryPart >= 0.0)
   {
      output << "+ " << z.mImaginaryPart << "i)";
   }
   else
   {
      // z.mImaginaryPart < 0.0
      // Replace + with minus sign
      output << "- " << -z.mImaginaryPart << "i)";
   }
}
//Code from Chapter06.tex line 779 save as ComplexNumber.cpp

//Methods called GetRealPart and GetImaginaryPart that allow us to
//access the corresponding private members.
double ComplexNumber::GetRealPart() const
{
    return mRealPart;
}
double ComplexNumber::GetImaginaryPart() const
{
    return mImaginaryPart;
}

//Friend functions RealPart and ImaginaryPart
//so one may either write z.GetImaginaryPart() or ImaginaryPart(z).
double RealPart(const ComplexNumber& z){
    return z.mRealPart;
}
double ImaginaryPart(const ComplexNumber& z){
    return z.mImaginaryPart;
}

//an overridden copy constructor
ComplexNumber::ComplexNumber(const ComplexNumber& z)
{
   mRealPart = z.mRealPart;
   mImaginaryPart = z.mImaginaryPart;
}

//a constructor that allows to specify a real number in complex form
//through a constructor that accepts 1 double prec. floating point as input
//sets real part of complex number to input and Im part to 0
ComplexNumber::ComplexNumber(double real)
{
    mRealPart = real;
    mImaginaryPart = 0.0;
}
//const method CalculateConjugate to get complex conjugate of input
ComplexNumber ComplexNumber::CalculateConjugate() const
{
   ComplexNumber z(mRealPart,-mImaginaryPart);
   return z;
}
//method SetConjugate which has a void return type
//and sets the complex number x + iy to its complex conjugate x - iy.
void ComplexNumber::SetConjugate()
{
    mImaginaryPart = -mImaginaryPart;
}

