#include <iostream>
#include <cmath>
#include <array>
/*Write code that dynamically allocates memory for three 2 × 2 matrices
of double precision floating point numbers, A, B, C,
and assigns values to the entries of A and B.
Let C = A + B.
Extend your code so that it calculates the entries of C,
and then prints the entries of C to screen.
Finally, de-allocate memory.
Again, check you have de-allocated memory correctly
by using a for loop as in the previous exercise. */

int main(){
for (int k = 0; k < 1e2; k++){
    int rows = 2, cols = 2;
    double** A;
    double** B;
    double** C;
    A = new double* [rows];
    B = new double* [rows];
    C = new double* [rows];
    for (int i = 0; i < rows; i++){
        A[i] = new double [cols];
        B[i] = new double [cols];
        C[i] = new double [cols];
    }
    for (int i = 0; i< 2; i++){
            for (int j = 0; j< 2; j++){
                A[i][j] = 1.0;
                B[i][j] = 3.0;
            }
    }
    for (int i = 0; i < 2; i++)
    {
        for(int j=0; j <2; j++){
        C[i][j] = A[i][j] + B[i][j];
        std::cout << C[i][j] << "\n";
        }
    }

    for (int i = 0; i < rows; i++)
    {
        delete[] A[i];
        delete[] B[i];
        delete[] C[i];
    }

    delete [] A;
    delete [] B;
    delete [] C;
    }
    return 0;

}
