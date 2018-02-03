#include <iostream>
#include <cmath>
double det(double **A);
int main(){
}
double det(double **A){
        double det_value;
        det0 = A[0][0]*( (A[1][1]*A[2][2]) - (A[1][2]*A[2][1]));
        det1 = A[0][1]*( (A[1][0]*A[2][2]) - (A[1][2]*A[2][0]));
        det2 = A[0][2]*( (A[1][0]*A[2][1]) - (A[1][1]*A[2][0]));
        det_value = det0 - det1 + det2;
        return det_value;
    }
