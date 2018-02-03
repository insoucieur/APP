#include "5_10.h"
#include <iostream>
int main(){
    int ARows = 3, ACols = 3;
    double** A = new double* [ARows];
    double* b = new double [ARows];
    double* u = new double [ARows];

    for (int i = 0; i< ARows; i++ )
    {
        A[i] = new double [ACols];
    }

    //assign values
    /*for (int i = 0; i< ARows; i++ ){
    for (int j = 0; j< ACols; j++ ){
        A[i][j] = ;
    }
    for (int i = 0; i< ARows; i++ ){
    b[i] = 9.0;
}*/
    A[0][0] = 1.0; A[0][1] = 1.0; A[0][2] = 1.0;
    A[1][0] = 2.0; A[1][1] = 1.0; A[1][2] = 3.0;
    A[2][0] = 3.0; A[2][1] = 1.0; A[2][2] = 6.0;
    b[0] = 4.0; b[1] = 7.0; b[2] = 2.0;

    guassian_elimination(A, b, u, 3);
    for (int i = 0; i< ARows; i++ ){
        std::cout << u[i] << " ";
}
for (int i = 0; i<ARows; i++)
    {
        delete[] A[i]; //delete array within matrix
    }

    delete[] A;
    delete[] b;
    delete[] u;

    //std::cout << b[0] << "\n";
    std::cout << "\n" << std::endl; }
