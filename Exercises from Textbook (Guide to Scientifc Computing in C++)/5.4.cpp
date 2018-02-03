#include <iostream>
#include <cmath>
/* function that can be used to calculate the mean and
standard deviation of an array of double precision
floating point numbers. */
void calculateMeanStd(int size, double* x, double& mean, double& stdev);
int main(){
    int n = 3;
    double* x = new double [n];
    for (int i = 0; i< 3; i++){
            x[i] = 3.0*double(i);
            std::cout << x[i] << "\n";
    }
    double mean, stdev;
    calculateMeanStd(n, x, mean, stdev);
    std::cout << "Mean is " << mean << "\n" << "Standard deviation is " << stdev << "\n";
    delete [] x;
    return 0;
}


void calculateMeanStd(int size, double* x, double& mean, double& stdev)
{

    double sum = 0.0;
    for (int i = 0; i< length; i++)
    {
        sum += a[i];}
    mean = sum/double(length);

    double stdev_sum = 0.0;
    for (int i = 0; i< length; i++)
    {
        stdev_sum += pow((a[i] - mean), 2);}

    if (length == 1){
        stdev = sqrt(stdev_sum/ double(length));
        }
    else{
    stdev = sqrt(stdev_sum/ double(length - 1));
    }
}
