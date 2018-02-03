#include <iostream>
#include <cmath>

int changeIntegerValue (int* p_i);
/* Write code that sends the address of an integer to a function
that changes the value of the integer. */
int main(){
    int i = 7;
    int* p_i = &i;
    i = changeIntegerValue(p_i);
    std::cout << i << "\n";
    return 0;
}

int changeIntegerValue(int* p_i)
{
    *p_i = 3;
    return *p_i;
}

