#include "5_3.h"
#include <iostream>

void swap_pointer(double* a, double* b){
    double c = *a;
    *a = *b;
    *b = c;
}
void swap_ref(double& a, double& b){
    double c = a;
    a = b;
    b = c;
}

