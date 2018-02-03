#include "6.1.hpp"
#include <cmath>
#include <iostream>

int main()
{
    ComplexNumber z1(0.0, 0.0);
    //std::cout << "z1 = " << z1 << "\n"; //does not
    std::cout << z1 << "\n"; // gives segmentation error when printed alone

    //printing both does not
}
