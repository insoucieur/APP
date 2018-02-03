#include <iostream>
using namespace std;
int main()
{
    /* 1. calculates sum of a collection of +ve integers inputted by user
    enter each integer followed by return key
    enter "-1" at the end of list of integers to be added
     */

     int sum = 0;
     int i = 0;
     while(i != -1)
     {
        cout << "Enter an integer \n and -1 at end of integers list \n and -2 to reset sum to zero \n";
        cin >> i;
        if (i != -1){
            sum += i;
        }
        else{
            cout << "Total sum is \n" << sum << "\n";
            break;
            }
     }
}
