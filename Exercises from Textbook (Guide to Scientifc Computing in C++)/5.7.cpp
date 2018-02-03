#include <iostream>
#include <cmath>
#include <cassert>
/* Calculate p-norm of a given vector, where p takes detault value 2
p-norm = [ sum_i (|x_i|**p ]** (1/p) */
double CalculateNorm (double* x, int vecSize, int p);
int main()
{
    int vecSize = 3;
    double* x = new double [vecSize];
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    int p = 2;
    double norm = CalculateNorm(x, vecSize, p);
    std::cout << norm << "\n";
}
double CalculateNorm (double* x, int vecSize, int p){
    double sum = 0.0;
    //Loop over eleements x_i of x, incremeting sum by |x_i|**p
    for (int i = 0; i < vecSize; i++)
    {
        double temp = fabs(x[i]);
        sum += pow(temp, p);

    }
    return pow(sum, 1.0/p);
}
