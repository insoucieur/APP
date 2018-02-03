#include "5_4.h"
#include <iostream>
#include <cmath>
double calc_std(double a[], int length){
    double mean;
    double sum = 0.0;
    for (int i = 0; i< length; i++)
    {
        sum += a[i];}
    mean = sum/double(length);
    double stdev;
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

    return stdev;
}

double calc_mean(double a[], int length){
    double mean;
    double sum = 0.0;
    for (int i = 0; i< length; i++)
    {
        sum += a[i];
}
    mean = sum/double(length);
    return mean;
}
