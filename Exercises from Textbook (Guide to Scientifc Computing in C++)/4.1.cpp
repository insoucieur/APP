#include <iostream>
#include <cmath>
/* integer i to take value 5 and declare pointer to integer p_j
store address of i in this pointer
multiply value of i by 5 using line that only uses pointer variable
declare another pointer to integer p_k
and use new keyword to allocate location in memory that this pointer stores
then store contents of variable i in this location */
int main(){
    int i = 5;
    int* p_j = &i; //p_j stores address of i
    int j = *p_j*5; //multiply i by 5 using pointer variable
    int* p_k;
    p_k = new int;
    *p_k = i;
    delete p_k;


}
