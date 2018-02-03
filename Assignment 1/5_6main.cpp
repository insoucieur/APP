#include "5_6.h"
#include <iostream>
#include <cassert>
int main(){
    int ARows = 10, ACols = 10;
    int BRows = 10, BCols = 10;
    double** A = new double* [ARows];
    double** res = new double* [BRows];
    double** B = new double* [BRows];
    double* a = new double [ARows];
    double* res_v = new double [ARows];
    double* v = new double [ARows];
    for (int i = 0; i< ARows; i++ )
    {
        A[i] = new double [ACols];
    }
    for (int i = 0; i< BRows; i++ )
    {
        B[i] = new double [BCols];
        res[i] = new double [BCols];
    }
    //assign values
   for (int i = 0; i< ARows; i++ ){
        for (int j = 0; j< ACols; j++ ){
            A[i][j] = 1.0;
        }
    }
    for (int i = 0; i< BRows; i++ ){
    for (int j = 0; j< BCols; j++ ){
        B[i][j] = 3.0;
    }
}
for (int i = 0; i< ARows; i++ ){
    a[i] = 1.0;
}

Multiply(res, A, B, 2, 3, 3, 2);
    for (int i = 0; i< 2; i++ ){
    for (int j = 0; j< 2; j++ ){
        std::cout << res[i][j] << " ";
    }
}
 std::cout << "\n";
Multiply(res_v, a, B, 3, 3, 2);
    for (int i = 0; i< 2; i++ ){
        std::cout << res_v[i] << " ";
}
std::cout << "\n";
Multiply(v, B, a, 2, 3, 3);
    for (int i = 0; i< 2; i++ ){
        std::cout << v[i] << " ";
}
std::cout << "\n";
double scalar = 3.0;
Multiply(res, scalar, A, 3, 3);
    for (int i = 0; i< 3; i++ ){
    for (int j = 0; j< 3; j++ ){
        std::cout << res[i][j] << " ";
    }
}
std::cout << "\n";
double scalar2 = 9.0;
Multiply(res, scalar2, A , 3, 3);
    for (int i = 0; i< 3; i++ ){
    for (int j = 0; j< 3; j++ ){
        std::cout << res[i][j] << " ";
    }
}

for (int i = 0; i<ARows; i++)
    {
        delete[] A[i]; //delete array within matrix
    }
for (int i = 0; i<BRows; i++)
    {
        delete[] B[i]; //delete array within matrix
        delete[] res[i];
    }
    delete[] A;
    delete[] res;
    delete[] B;
    delete[] a;
    delete[] res_v;
    delete[] v;

    std::cout << res[0][0] << "\n";
    std::cout << "\n" << std::endl; }
