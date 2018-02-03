#include <iostream>
#include <cmath>
/* Assign values to 2 integer variables.
Swap values stored by these variables using only pointers to integers. */
int main(){
    int a = 3;
    int b = 7;

    int* p_a = &b;
    int* p_b = &a;
    //*p_a = *p_b;
    //*p_b = *p_a;
    //*p_a is contents of memory
    std::cout << *p_a << "\n";
    std::cout << *p_b;

}
