#include <iostream>
#include <cmath>
using namespace std;
int main (int argc, char* argv[])
{
    double u[3] = {1.0, 2.0, 3.0};
    double v[3] = {6.0, 5.0, 4.0};
    double A[3][3] = {{1.0, 5.0, 0.0},
                    {7.0, 1.0, 2.0},
                    {0.0, 0.0, 1.0}};
    double B[3][3] = {{-2.0, 0.0, 1.0},
                {1.0, 0.0, 0.0},
                {4.0, 1.0, 0.0}};
    double w[3];
    for (int i=0; i < 3; i++)
    {
        w[i] = u[i] - 3.0*v[i];

    }
    /* x = u - v
        y = Au
        z = Au -v
        C = 4A -3 B
        D = AB
        [rows][columns]
     */
    double x[3];
    for (int i = 0; i<3; i++)
    {
        x[i] = u[i] - v[i];
    }
    double y[3];
    for (int i = 0; i< 3; i++)
    {
        for(int j=0; j <3; j++){
        y[i] += A[i][j]*u[j];
        }
    }
    double z[3];
    for (int i = 0; i<3; i++)
    {
        z[i] += y[i] - v[i];
    }
    double C[3][3];
    for (int i = 0; i < 3; i++)
    {
        for(int j=0; j <3; j++){
        C[i][j] = 4*A[i][j] - 3*B[i][j];
        }
    }
    double D[3][3];
    for (int a = 0; a < 3; a++)
    {
        for(int b = 0; b < 3; b++){
            for (int i = 0; i < 3; i++){
            D[a][b] += A[a][i]* B[i][b];
            }
        }
    }
    cout << D[0][0] << "\n";
    cout << D[2][1] << "\n";
    cout << C[0][0] << "\n";
    cout << C[2][1] << "\n";
    return 0;
}
