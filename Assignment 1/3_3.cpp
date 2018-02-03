#include "3_3.h"
#include <iostream>
#include <cmath>
#include <fstream>
void implicit_Euler(int n){
     //use N to calculate h
     double const h = double (1)/ double (n);//step size
     //x = nh for n = 0, 1, 2,.. N-1
     //xy.dat has 2 columns, calculated values of x and calculated values of y
     //(y_n - y_n-1)/h = -y_n with y_0 = 1
    std::ofstream write_output ("xy.dat");
    double y_prev = 1.0;
    for (int i = 0; i< n; i++){
        double x_i = double(i)*h;
        double y_i = y_prev;
        write_output << x_i << "," << y_i << "\n";
        y_i = y_prev / (1.0 + h);
        y_prev = y_i;
    }
    write_output.close();

}


