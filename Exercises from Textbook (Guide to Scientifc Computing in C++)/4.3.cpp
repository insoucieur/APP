#include <iostream>
#include <cmath>
/*Write code that allocates memory dynamically to two vectors
of double precision floating point numbers of length 3,
assigns values to each of the entries,
and then de-allocates the memory before the code terminates.
Extend this code so that it calculates the scalar (dot) product of these vectors
and prints it to screen before the memory is de-allocated.
Put the allocation of memory, calculation and de-allocation of
memory inside a for loop that runs 1,000,000,000 times:
if the memory is not de-allocated properly your code will
use all available resources and your computer may struggle. */

int main(){
    for (int j = 0; j < 1e3; j++){
        double* x;
        double* y;
        x = new double [3];
        y = new double [3];

        for (int i = 0; i< 3; i++){
            x[i] = double(i);
            y[i] = 3.0*x[i];
        }
        delete[] x;
        delete[] y;
    }
    return 0;
}
