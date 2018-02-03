#include "5_6.h"
#include <iostream>
#include <cmath>
void Multiply(double **res, double **A, double **B, int ARows, int ACols, int BRows, int BCols){
    if (ACols == BRows){
        for (int a = 0; a < ARows; a++){
        for(int b = 0; b < BCols; b++){
            for (int i = 0; i < ACols; i++){
            res[a][b] += A[a][i]* B[i][b];
            }
            }
        }
    }
}

void Multiply(double *res, double *A, double **B, int ACols, int BRows, int BCols){
    if (ACols == BRows){
        for(int r = 0; r < BRows; r++)
        {
            for (int c = 0; c < BCols; c++){
            res[c] += A[r]*B[r][c];
            }}
    }
}
void Multiply(double *res, double **A, double *B, int ARows, int ACols, int BRows){
    if (BRows == ACols){
    for(int r = 0; r < ARows; r++)
    {
        for (int c = 0; c < ACols; c++){
        res[r] += A[r][c]*B[c];
        }}
        }
}
void Multiply(double **res, double scalar, double **B, int BRows, int BCols){
    for(int r = 0; r < BRows; r++)
    {
        for (int c = 0; c < BCols; c++){
        res[r][c] = scalar*B[r][c];
        }
        }
}
void Multiply(double **res, double **B, double scalar, int BRows, int BCols){
    for(int r = 0; r < BRows; r++)
    {
        for (int c = 0; c < BCols; c++){
        res[r][c] = B[r][c]*scalar;
        }
        }
}
