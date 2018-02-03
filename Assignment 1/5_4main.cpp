#include "5_4.h"
#include <iostream>
int main(){
    int n = 3;
    double* x = new double [n];
    for (int i = 0; i< n; i++){
            x[i] = 3.0*double(i);
            }
    std::cout << calc_std(x, n) << "\n";
    std::cout << calc_mean(x, n)<< std::endl;

    delete[]x; }
