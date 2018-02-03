#include <iostream>
#include <cmath>

void printIntegerValue (int* p_i);
/* Write code that sends the address of an integer to a function
that prints out the value of the integer. */
int main(){
    int i = 7;
    int* p_i = &i;
    printIntegerValue(p_i);

    return 0;
}

void printIntegerValue(int* p_i)
{
    std::cout << *p_i << "\n";
}

