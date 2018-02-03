#include "5_3.h"
#include <iostream>
int main(){
    double a = 7.0;
    double b = 3.0;
    swap_pointer(&a, &b);
    std::cout << a << " " << b << "\n";
    swap_ref(a, b);
    std::cout << a << " " << b << std::endl; }


