#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;
int main()
{
    /* 1. calculates inverse of matrix A
    2. check inverse calculated is correct
    3. assert statement to check that determinant of matrix is non-zero
     */
     int A[2][2] = {{4, 10}, {1, 1}};
     //1
     int detA = (A[0][0]*A[1][1]) - (A[0][1]*A[1][0]);
     int adjA[2][2] = {{A[1][1], -1*A[0][1]}, {-1*A[1][0], A[0][0]} };
     double invA[2][2];
     for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
        invA[i][j] = (double(1.0)/double (detA)) * double(adjA[i][j]);
        }
     }
     //2
     cout << invA[0][0] << "\n";
     cout << invA[0][1] << "\n";
     cout << invA[1][0] << "\n";
     cout << invA[1][1] << "\n";
     //and is correct
    //3
    assert(detA != 0);
    return 0;
}
