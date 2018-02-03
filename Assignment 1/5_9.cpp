#include "5_9.h"
#include <iostream>
#include <cmath>

void solve3by3 ( double **A, double *b , double *u){
    //make a copy of the matrix A and vector b
    double mat[3][3];
    double vec[3];
    for (int i = 0; i< 3; i++ ){
            vec[i] = b[i];
        for (int j = 0; j< 3; j++ ){
            mat[i][j] = A[i][j];
        }
    }
    //calculate determinant
    double det, det0, det1, det2;
    det0 = mat[0][0]*( (mat[1][1]*mat[2][2]) - (mat[1][2]*mat[2][1]));
    det1 = mat[0][1]*( (mat[1][0]*mat[2][2]) - (mat[1][2]*mat[2][0]));
    det2 = mat[0][2]*( (mat[1][0]*mat[2][1]) - (mat[1][1]*mat[2][0]));
    det = det0 - det1 + det2;

    //calculate cofactor
    double cofactor[3][3];
    cofactor[0][0] = ((mat[1][1]*mat[2][2]) - (mat[1][2]*mat[2][1]));
    cofactor[1][0] = ((mat[1][2]*mat[2][0]) - (mat[1][0]*mat[2][2]));
    cofactor[2][0] = ((mat[1][0]*mat[2][1]) - (mat[1][1]*mat[2][0]));

    cofactor[0][1] = ((mat[0][2]*mat[2][1]) - (mat[0][1]*mat[2][2]));
    cofactor[1][1] = ((mat[0][0]*mat[2][2]) - (mat[0][2]*mat[2][0]));
    cofactor[2][1] = ((mat[0][1]*mat[2][0]) - (mat[0][0]*mat[2][1]));

    cofactor[0][2] = ((mat[0][1]*mat[1][2]) - (mat[0][2]*mat[1][1]));
    cofactor[1][2] = ((mat[0][2]*mat[1][0]) - (mat[0][0]*mat[1][2]));
    cofactor[2][2] = ((mat[0][0]*mat[1][1]) - (mat[0][1]*mat[1][0]));

    double inv_mat[3][3];
    for (int i=0; i < 3; i++){
        for (int j=0; j < 3; j++){
        inv_mat[i][j] = cofactor[i][j]/det;
    }}

    for(int r = 0; r < 3; r++)
    {
        for (int c = 0; c < 3; c++){
        u[r] += inv_mat[r][c]*vec[c];
        }}


}
