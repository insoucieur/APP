#include <iostream>
#include <cmath>
using namespace std;
int main()
{ /* 1. || = OR
    (x > y) || (x < 5.0) => if x is greater than y or x is less than 5
    then z = 4.0
    otherwise z = 2.0

    2. a) z = 4.0
    b) z = 2.0
    c) z = 4.0

    3. replace code so x > y is replaced by x >= y */
    double x = 7.0, y = 10.0, z;
    if ((x >= y) || (x < 5.0))
    {
        z = 4.0;
    }
    else
    {
        z = 2.0;
    }
    cout << z;
}

