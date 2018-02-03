#include <iostream>
#include <cmath>


/* Write a function that swaps the values of two double precision floating
point numbers, so that these changes are visible in the code that has
called this function.
1. Write this function using pointers.
2. Write this function using references. */
void swapValues(double *p_a, double *p_b);
void swapUsingRef(double& a, double &b);
int main(){
    double a = 7.0;
    double b = 3.0;
    //swapValues(&a, &b);
    swapUsingRef(a, b);
    std::cout << a << "\n";
    std::cout << b << "\n";
    return 0;
}

void swapValues(double *p_a,  double *p_b)
{
    double c = *p_a; //has value of a
    *p_a = *p_b; //b takes value of a
    *p_b = c; //b takes value of a

}

void swapUsingRef(double& a, double& b)
{
    double c = a;
    a = b;
    b = c;
}

