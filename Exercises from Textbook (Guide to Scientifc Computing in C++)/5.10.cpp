#include "5_10.h"
#include <iostream>
#include <cmath>

void guassian_elimination(double **A, double *b, double *u, int n){
    //make a copy of the matrix A and vector b
    double X[n][n];
    double y[n];
    for (int i = 0; i< n; i++ ){
            y[i] = b[i];
        for (int j = 0; j< n; j++ ){
            X[i][j] = A[i][j];
        }
    }
    //find row n with largest absolute value of a_nk (so such through all
    //elements of each k column in the kth step for n = k, k+1, ..., N-1,
    // then interchange that row with row k using the P matrix
    //until abs(a_kk) is not 0
    double lambda_r;
    for (int k = 0; k < n-1; k++){
        // temp value to initialise loop to find kmax
        double mmax_X = fabs(X[k][k]);
        int mmax = k;
        for (int m = k+1; m < n; m++){
            if ((fabs(X[m][k])) > mmax_X){
                //loop through each m and replace the last biggest value
                //with the current one until the largest is found
                mmax_X = fabs(X[m][k]);
                mmax = m;
            }}
            //std::cout << "m_max " << mmax << "\n";
        //swap row with max with current k row, need to swap rows of vector as well!
        for (int c = 0; c < n; c++){
            double X_rowmax = X[mmax][c];
            X[mmax][c] = X[k][c];
            X[k][c] = X_rowmax;
            //std::cout << X[k][c] << "\n";
        }
        double y_rowmax = y[mmax];
        y[mmax] = y[k];
        y[k] = y_rowmax;
        //std::cout << y[k] << "\n";
    //then need to zero the rows for each column k below the diagonal to
    //and form an upper triangular matrix after interchanging rows
    for (int r = k+1; r < n; r++){
        lambda_r = X[r][k]/X[k][k];
        //std::cout << "vector before " << y[r] << "\n";
        //std::cout << lambda_r << y[k] << "\n";
        y[r] = y[r] - (lambda_r*y[k]);
        //std::cout << "vector after " << y[r] << "\n";
        //std::cout << "lambda " << lambda_r << "\n";
        for (int c = 0; c < n; c++){
        //std::cout << "before " << X[r][c] << "\n";
        X[r][c] = X[r][c] - (lambda_r*X[k][c]);
        //std::cout << "after " << X[r][c] << "\n";
        }}
        }
    //then solve the upper triangular system from the last eqn in the system
    // so need to iterate backwards from N-1 to 0
    for (int k = n-1; k>=0; k--){
        double sum = 0.0;
        for (int i = k+1; i< n; i++){
                sum += X[k][i]*u[i];
        }
        u[k] = (y[k] - sum)/X[k][k];
    }
        }
