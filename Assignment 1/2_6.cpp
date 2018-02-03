#include "2_6.h"
#include <iostream>
#include <cmath>
double newton_Raphson(double initialGuess, double epsilon){
    double x_prev = initialGuess;
    double x_next = x_prev - ((exp(x_prev) + pow(x_prev,3) -5)/(exp(x_prev) + 3*pow(x_prev,2)));
    while (fabs(x_next - x_prev) > epsilon){
        x_prev = x_next;
        x_next = x_prev - ((exp(x_prev) + pow(x_prev,3) -5)/(exp(x_prev) + 3*pow(x_prev,2)));
    }
    return x_next;
}

