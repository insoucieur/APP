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
        //for a diagonal matrix
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
