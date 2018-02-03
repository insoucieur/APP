#include <iostream>
#include <fstream>
int main (int arg, char* argv[])
{
    double x[4] = {0.0, 1.0, 1.0, 0.0};
    double y[4] = {0.0, 0.0, 1.0, 1.0};


    /* extend code to print arrays x and y called x_and_y.dat
    so data file has 4 elements of x on top line and 4 elements of y on next line
    2. output sream flushed immediately after each line of file is written
    3. precision is set to 10 s.f. and output is in scientific notation & plus signs for +ve numbers
    */
    std::ofstream write_output ("x_and_y.dat");
    write_output.setf(std::ios::scientific); //output in scientific format
    write_output.setf(std::ios::showpos); //always show + or - sign
    write_output.precision(10); //10 sig figs

    for (int i = 0; i< 4; i++)
    {
        write_output << x[i] << " ";
    }
    write_output << "\n";
    write_output.flush(); //2
    for (int i = 0; i< 4; i++){
    write_output <<  y[i] << " ";
    }
    write_output.flush(); //2
    write_output.close();
    return 0;
}

