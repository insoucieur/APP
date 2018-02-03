/* Additional general notes */
-> void => a way to calculate return multiple values since double functions can only return one
-> new always followed by delete
-> if function e.g. Student(std::string name) contains parameter that has to be assigned to a variable of the same name use:
	this->name = name;

double b; //declare variable; allocate memory for a double precision floating point
b = 3.0; //set value

int a = 7, c = 0; //can also initialise when defining variable typ(variables may not be initialised to 0 when memory is allocated)

const double d = 10.0; //variable intended to be constant - guaranteed to be unchanged

#include <iostream>
std::cout << a << “\n”; //to print a
std::cout.flush(); //make sure output a is printed before executing other statements

//keyboard input:
int pin, account_number;
std::cout << “Enter account number \n”;
std::cout << “Enter PIN, then hit return \n”;
std::cin >> account_number >> pin;

//for strings
std::string name;
std::out << “Enter name then hit return \n”;
std::getline(std::cin, name);
std::cout << “Your name is “ << name << \n”;

/* Basic maths */
#include <cmath>

double x, y;
sqrt(x) //square root
exp(x) //exponential fn
pow(x, y) //rise x to power of y
M_PI //pi

a+=b; //a=a+b
a-=b; //a=a-b
a*=b; //a=a*b
a/=b; //a=a/b
a%=b; //a=a%b
a++; //a=a+1
a—-; //a=a-1

int i = 5, j =2;
double k;
k = ( (double)(i)/ (double)(j) ) //explicit type conversion - convert integers to double

/* ——- */

//indices of an n sized array start at 0 and end at n-1
//initialising arrays:
int array1[2]; //vector of integers of length 2
double array2[2][3]; //matrix of double of size 2 by 3
//then set values of arrays:
array1[0] = 1;
array1[1] = 3;

//or initialise arrays when declared:
double array3[3] = {5.0, 1.0, 2.0};
int array4[2][2] = { {1,6,-4}, {2,2,2} };

array1[0]++ //increments value of this entry by 1

//will not be accepted by the code:
int array6[3];
array6 = {0, 1, 2}; //cannot assign values like this after it’s been initialised
int array5[3] = {0,1,2} //but can be declared this way

char letter;
letter = ‘a’; //single quotation marks
std::string city;
city = “Oxford”; //double quotation marks
std::cout << “String length = “ << city.length() << “\n”;
std::cout << “Third character = “ >> city.at(2) << “\n”; //can also do city[2]
std::cout << city << “\n”; //prints the string in city, can also do city.c_str()


/* if, else, else if, while, for statements */
//if, else, else if:
double x = -2.0;
if ( (x < 0.0) && (x%2 == 0) ) //&& = AND, || = or ! => NOT
{
	x = 0.0;
}
else if (x > 1.0)
{
	x = 1.0;
}
else
{
	x = 2.0;
}

//while e.g. function for Newton_Raphson method:
double newton_Raphson(double initialGuess, double epsilon){
    double x_prev = initialGuess; //x_0
    double x_next = x_prev - ((exp(x_prev) + pow(x_prev,3) -5)/(exp(x_prev) + 3*pow(x_prev,2))); //get y_1 using x_0
    while (fabs(x_next - x_prev) > epsilon){ //while this is true, do iteration to get y_n from x_n-1
        x_prev = x_next; //x_n-1 is now the last calculated value x_n
        x_next = x_prev - ((exp(x_prev) + pow(x_prev,3) -5)/(exp(x_prev) + 3*pow(x_prev,2))); //calculate x_n until x_n ~ x_prev accurate to precision specified by epsilon
    }
    return x_next;
}

/* ——- */

/* ASSIGNMENT 1*/

//writing output to file
#include <iostream>
#include <cmath>
#include <fstream>
void implicit_Euler(int n){
     //use N to calculate h
     double const h = double (1)/ double (n);//step size
     //x = nh for n = 0, 1, 2,.. N-1
     //xy.dat has 2 columns, calculated values of x and calculated values of y
     //(y_n - y_n-1)/h = -y_n with y_0 = 1
    std::ofstream write_output ("xy.dat");
    double y_prev = 1.0;
    for (int i = 0; i< n; i++){
        double x_i = double(i)*h; //this
        double y_i = y_prev; //and this defined so first value is written in file
        write_output << x_i << "," << y_i << "\n"; //writes value of x and y after each calculation
        y_i = y_prev / (1.0 + h); //calculate next value
        y_prev = y_i; //calculated value is now used in the next one until n-1 is reached
    }
    write_output.close();

}

//using pointers and references
void swap_pointer(double* a, double* b){
    double c = *a; //should be double* c = *a
    *a = *b;
    *b = c; //and *b = *c
}
void swap_ref(double& a, double& b){
    double c = a;
    a = b;
    b = c;
}
//calculate mean and std dev
#include "5_4.h"
#include <iostream>
#include <cmath>
//a[] is a pointer to an array
double calc_std(double a[], int length){ // array a of unspecified length that will be specified by int length variable
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
        stdev_sum += pow((a[i] - mean), 2);} //the sum part in square root
    //distinguish between if it is length 1 or not so division can be carried out and then square root
    if (length == 1){
        stdev = sqrt(stdev_sum/ double(length)); //avoid division by zero when length is 1
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

//multiplying matrices and vectors using pointers
//function to multiply a matrix and a vector
void Multiply(double *res, double **A, double *B, int ARows, int ACols, int BRows){ //by writing fn implementation before use i.e. int main()
    //can skip function prototype e.g. void Multiply(double *res, double **A, double *B, int ARows, int ACols, int BRows); before int main() then specifying implementation after
    if (BRows == ACols){
    for(int r = 0; r < ARows; r++)
    {
        for (int c = 0; c < ACols; c++){
        res[r] += A[r][c]*B[c];
        }}
        }
}
int main(){
    int ARows = 10, ACols = 10;
    int BRows = 10, BCols = 10;
    double** A = new double* [ARows];
    double** res = new double* [BRows];
    double** B = new double* [BRows];
    double* a = new double [ARows];
    double* res_v = new double [ARows];
    double* v = new double [ARows];
    for (int i = 0; i< ARows; i++ )
    {
        A[i] = new double [ACols];
    }
    for (int i = 0; i< BRows; i++ )
    {
        B[i] = new double [BCols];
        res[i] = new double [BCols];
    }
    //assign values
   for (int i = 0; i< ARows; i++ ){
        for (int j = 0; j< ACols; j++ ){
            A[i][j] = 1.0;
        }
    }
    for (int i = 0; i< BRows; i++ ){
    for (int j = 0; j< BCols; j++ ){
        B[i][j] = 3.0;
    }
}
for (int i = 0; i< ARows; i++ ){
    a[i] = 1.0;
}


std::cout << "\n";
Multiply(v, B, a, 2, 3, 3);
    for (int i = 0; i< 2; i++ ){
        std::cout << v[i] << " ";
}

//free allocated memory
for (int i = 0; i<ARows; i++)
    {
        delete[] A[i]; //delete array within matrix
    }
for (int i = 0; i<BRows; i++)
    {
        delete[] B[i]; //delete array within matrix
        delete[] res[i];
    }
    delete[] A;
    delete[] res;
    delete[] B;
    delete[] a;
    delete[] res_v;
    delete[] v;

    std::cout << "\n" << std::endl; }

//to solve a 3x3 linear system
void solve3by3 ( double **A, double *b , double *u){
    //make a copy of the matrix A and vector b
    double mat[3][3];
    double vec[3];
    for (int i = 0; i< 3; i++ ){
            vec[i] = b[i];
        for (int j = 0; j< 3; j++ ){
            mat[i][j] = A[i][j];
        }
    }
    //calculate determinant
    double det, det0, det1, det2;
    det0 = mat[0][0]*( (mat[1][1]*mat[2][2]) - (mat[1][2]*mat[2][1]));
    det1 = mat[0][1]*( (mat[1][0]*mat[2][2]) - (mat[1][2]*mat[2][0]));
    det2 = mat[0][2]*( (mat[1][0]*mat[2][1]) - (mat[1][1]*mat[2][0]));
    det = det0 - det1 + det2;

    //calculate cofactor
    double cofactor[3][3];
    cofactor[0][0] = ((mat[1][1]*mat[2][2]) - (mat[1][2]*mat[2][1]));
    cofactor[1][0] = ((mat[1][2]*mat[2][0]) - (mat[1][0]*mat[2][2]));
    cofactor[2][0] = ((mat[1][0]*mat[2][1]) - (mat[1][1]*mat[2][0]));

    cofactor[0][1] = ((mat[0][2]*mat[2][1]) - (mat[0][1]*mat[2][2]));
    cofactor[1][1] = ((mat[0][0]*mat[2][2]) - (mat[0][2]*mat[2][0]));
    cofactor[2][1] = ((mat[0][1]*mat[2][0]) - (mat[0][0]*mat[2][1]));

    cofactor[0][2] = ((mat[0][1]*mat[1][2]) - (mat[0][2]*mat[1][1]));
    cofactor[1][2] = ((mat[0][2]*mat[1][0]) - (mat[0][0]*mat[1][2]));
    cofactor[2][2] = ((mat[0][0]*mat[1][1]) - (mat[0][1]*mat[1][0]));

    //calculate inverse
    double inv_mat[3][3];
    for (int i=0; i < 3; i++){
        for (int j=0; j < 3; j++){
        inv_mat[i][j] = cofactor[i][j]/det;
    }}

    //use inverse to solve system
    for(int r = 0; r < 3; r++)
    {
        for (int c = 0; c < 3; c++){
        u[r] += inv_mat[r][c]*vec[c];
        }}

}

//Gaussian elimination with pivoting
void guassian_elimination(double **A, double *b, double *u, int n){
    //make a copy of the matrix A (n by n) and vector b (length n)
    double X[n][n];
    double y[n];
    for (int i = 0; i< n; i++ ){
            y[i] = b[i];
        for (int j = 0; j< n; j++ ){
            X[i][j] = A[i][j];
        }
    }
    //find row n with largest absolute value of a_nk (so sort through all
    //elements of each k column in the kth step for n = k, k+1, ..., N-1,
    // then interchange that row with row k using the P matrix
    //until abs(a_kk) is not 0
    double lambda_r;
    for (int k = 0; k < n-1; k++){
        // temp value to initialise loop to find kmax
        double mmax_X = fabs(X[k][k]);
        int mmax = k;
        for (int m = k+1; m < n; m++){
            if ((fabs(X[m][k])) > mmax_X){
                //loop through each m and replace the last biggest value
                //with the current one until the largest is found
                mmax_X = fabs(X[m][k]);
                mmax = m;
            }}
            //std::cout << "m_max " << mmax << "\n";
        //swap row with max with current k row, need to swap rows of vector as well!
        for (int c = 0; c < n; c++){
            double X_rowmax = X[mmax][c];
            X[mmax][c] = X[k][c];
            X[k][c] = X_rowmax;
            //std::cout << X[k][c] << "\n";
        }
        double y_rowmax = y[mmax];
        y[mmax] = y[k];
        y[k] = y_rowmax;
        //std::cout << y[k] << "\n";
    //then need to zero the rows for each column k below the diagonal to
    //and form an upper triangular matrix after interchanging rows
    for (int r = k+1; r < n; r++){
        lambda_r = X[r][k]/X[k][k];
        //std::cout << "vector before " << y[r] << "\n";
        //std::cout << lambda_r << y[k] << "\n";
        y[r] = y[r] - (lambda_r*y[k]);
        //std::cout << "vector after " << y[r] << "\n";
        //std::cout << "lambda " << lambda_r << "\n";
        for (int c = 0; c < n; c++){
        //std::cout << "before " << X[r][c] << "\n";
        X[r][c] = X[r][c] - (lambda_r*X[k][c]);
        //std::cout << "after " << X[r][c] << "\n";
        }}
        }
    //then solve the upper triangular system from the last eqn in the system
    // so need to iterate backwards from N-1 to 0
    for (int k = n-1; k>=0; k--){
        double sum = 0.0;
        for (int i = k+1; i< n; i++){
                sum += X[k][i]*u[i];
        }
        u[k] = (y[k] - sum)/X[k][k];
    }
        }
//to test it
#include <iostream>
int main(){
    int ARows = 3, ACols = 3;
    double** A = new double* [ARows];
    double* b = new double [ARows];
    double* u = new double [ARows];

    for (int i = 0; i< ARows; i++ )
    {
        A[i] = new double [ACols];
    }

    //assign values
    /*for (int i = 0; i< ARows; i++ ){
    for (int j = 0; j< ACols; j++ ){
        A[i][j] = ;
    }
    for (int i = 0; i< ARows; i++ ){
    b[i] = 9.0;
}*/
    A[0][0] = 1.0; A[0][1] = 1.0; A[0][2] = 1.0;
    A[1][0] = 2.0; A[1][1] = 1.0; A[1][2] = 3.0;
    A[2][0] = 3.0; A[2][1] = 1.0; A[2][2] = 6.0;
    b[0] = 4.0; b[1] = 7.0; b[2] = 2.0;

    guassian_elimination(A, b, u, 3);
    for (int i = 0; i< ARows; i++ ){
        std::cout << u[i] << " ";
}
for (int i = 0; i<ARows; i++)
    {
        delete[] A[i]; //delete array within matrix
    }

    delete[] A;
    delete[] b;
    delete[] u;

    //std::cout << b[0] << "\n";
    std::cout << "\n" << std::endl; }

/* Header files and constructors */

//GENERAL CASE:
//define header file = General.hpp
#ifndef GENERALHEADERDEF
#define GENERALHEADERDEF

#include <iostream>

class General
{
    private:
        double member1;
        double member2;
    public:
        General(); //default constructor
	//General::General(){ //set private member variables to 0.0 here}
        General(double x, double y); //constructor that sets values
	//General::General(double x, double y){//set private member variables using x and y here}
        double CalculateModulus() const; //method that returns a value
	//double General::CalculateModulus() const{}
	//initialise the function variable first then the class then the function
	General NegativeVersion() const; //method that returns a value of the same type (form) as the class General
	//General General NegativeVersion() const{}
        void SetNegative(); //method to set values without returning anything
	//void General::SetNegative(){}
        General(const General& z); //copy constructor
	//General::General(const General& z){member1 = z.member1; etc.}
	General(double onevalue); //constructor that specifies one value in same form as General
	//General::General(double onevalue){}

	//operators (see ComplexNumbers for how to define in .cpp):
        General& operator=(const General& z); //overload assignment (=) operator
        General operator-() const; //overload unary - operator
        General operator+(const General& z) const; //overload unary + operator
        General operator-(const General& z) const;
        friend std::ostream& operator<<(std::ostream& output,const General& z);

        //ways to access private members:
        double GetMember1() const; //method way to get private member variables
	//double General GetMember1() const{return member1;}
	//when using in main, use GetMember1(z)

        friend double Member2(const General& z); //friend function to get private member variables
	//double Member2(const General& z){return z.Member2;}
	//when using in main, have to use z.Member2

};

/* ASSIGNMENT 2*/

#ifndef COMPLEXNUMBERHEADERDEF
#define COMPLEXNUMBERHEADERDEF

#include <iostream>

class ComplexNumber
{
    private:
        double mRealPart;
        double mImaginaryPart;
    public:
        ComplexNumber();
        ComplexNumber(double x, double y);
        double CalculateModulus() const;
        double CalculateArgument() const;
        ComplexNumber CalculatePower(double n) const;
        ComplexNumber& operator=(const ComplexNumber& z);
        ComplexNumber operator-() const;
        ComplexNumber operator+(const ComplexNumber& z) const;
        ComplexNumber operator-(const ComplexNumber& z) const;
        friend std::ostream& operator<<(std::ostream& output,const ComplexNumber& z);

        //exercise prototypes
        double GetRealPart() const;
        double GetImaginaryPart() const;
        friend double RealPart(const ComplexNumber& z);
        friend double ImaginaryPart(const ComplexNumber& z);
        ComplexNumber(const ComplexNumber& z);
        ComplexNumber(double real);
        ComplexNumber CalculateConjugate() const;
        void SetConjugate();

        //not mandatory, but useful for exercise 6.1.7
        ComplexNumber operator*(const ComplexNumber& z) const;

};

#endif

//define functions, constructors in corresponding .cpp file
#include "ComplexNumber.hpp"
#include <cmath>
#include <iostream>

// Override default constructor
// Set real and imaginary parts to zero
ComplexNumber::ComplexNumber()
{
   mRealPart = 0.0;
   mImaginaryPart = 0.0;
}

// Constructor that sets complex number z=x+iy
ComplexNumber::ComplexNumber(double x, double y)
{
   mRealPart = x;
   mImaginaryPart = y;
}

// Method for computing the modulus of a
// complex number
double ComplexNumber::CalculateModulus() const
{
   return sqrt(mRealPart*mRealPart+
               mImaginaryPart*mImaginaryPart);
}

// Method for computing the argument of a
// complex number
double ComplexNumber::CalculateArgument() const
{
   return atan2(mImaginaryPart, mRealPart);
}

// Method for raising complex number to the power n
// using De Moivre's theorem - first complex
// number must be converted to polar form
ComplexNumber ComplexNumber::CalculatePower(double n) const
{
   double modulus = CalculateModulus();
   double argument = CalculateArgument();
   double mod_of_result = pow(modulus, n);
   double arg_of_result = argument*n;
   double real_part = mod_of_result*cos(arg_of_result);
   double imag_part = mod_of_result*sin(arg_of_result);
   ComplexNumber z(real_part, imag_part);
   return z;
}

// Overloading the = (assignment) operator
ComplexNumber& ComplexNumber::
               operator=(const ComplexNumber& z)
{
   mRealPart = z.mRealPart;
   mImaginaryPart = z.mImaginaryPart;
   return *this;
}

// Overloading the unary - operator
ComplexNumber ComplexNumber::operator-() const
{
   ComplexNumber w;
   w.mRealPart = -mRealPart;
   w.mImaginaryPart = -mImaginaryPart;
   return w;
}

// Overloading the binary + operator
ComplexNumber ComplexNumber::
              operator+(const ComplexNumber& z) const
{
   ComplexNumber w;
   w.mRealPart = mRealPart + z.mRealPart;
   w.mImaginaryPart = mImaginaryPart + z.mImaginaryPart;
   return w;
}

// Overloading the binary - operator
ComplexNumber ComplexNumber::
              operator-(const ComplexNumber& z) const
{
   ComplexNumber w;
   w.mRealPart = mRealPart - z.mRealPart;
   w.mImaginaryPart = mImaginaryPart - z.mImaginaryPart;
   return w;
}

// Overloading the insertion << operator
std::ostream& operator<<(std::ostream& output,
                         const ComplexNumber& z)
{
   // Format as "(a + bi)" or as "(a - bi)"
   output << "(" << z.mRealPart << " ";
   if (z.mImaginaryPart >= 0.0)
   {
      output << "+ " << z.mImaginaryPart << "i)";
   }
   else
   {
      // z.mImaginaryPart < 0.0
      // Replace + with minus sign
      output << "- " << -z.mImaginaryPart << "i)";

	return output; //NEEDED since std::ostream is kind of a type, so need to return something since it is not a void function
   }
}
//Code from Chapter06.tex line 779 save as ComplexNumber.cpp

//Methods called GetRealPart and GetImaginaryPart that allow us to
//access the corresponding private members.
double ComplexNumber::GetRealPart() const
{
    return mRealPart;
}
double ComplexNumber::GetImaginaryPart() const
{
    return mImaginaryPart;
}

//Friend functions RealPart and ImaginaryPart
//so one may either write z.GetImaginaryPart() or ImaginaryPart(z).
double RealPart(const ComplexNumber& z){
    return z.mRealPart;
}
double ImaginaryPart(const ComplexNumber& z){
    return z.mImaginaryPart;
}

//an overridden copy constructor
ComplexNumber::ComplexNumber(const ComplexNumber& z)
{
   mRealPart = z.mRealPart;
   mImaginaryPart = z.mImaginaryPart;
}

//a constructor that allows to specify a real number in complex form
//through a constructor that accepts 1 double prec. floating point as input
//sets real part of complex number to input and Im part to 0
ComplexNumber::ComplexNumber(double real)
{
    mRealPart = real;
    mImaginaryPart = 0.0;
}
//const method CalculateConjugate to get complex conjugate of input
ComplexNumber ComplexNumber::CalculateConjugate() const //ComplexNumber because return type is ComplexNumber
{
   ComplexNumber z(mRealPart,-mImaginaryPart);
   return z;
}
//method SetConjugate which has a void return type
//and sets the complex number x + iy to its complex conjugate x - iy.
void ComplexNumber::SetConjugate()
{
    mImaginaryPart = -mImaginaryPart;
}

//using ComplexNumber to calculate exponential of a complex diagonal matrix
//header file:
#ifndef _CALCULATEEXPONENTIAL_
#define _CALCULATEEXPONENTIAL_

#include "ComplexNumber.hpp"

void CalculateExponential(ComplexNumber **A, int nMax, ComplexNumber **res);


//non mandatory
void printMatrix(ComplexNumber **A, int rows, int cols);

#endif

//cpp file:
#include "CalculateExponential.hpp"
#include <cmath>
#include <iostream>
#include "ComplexNumber.hpp"
void CalculateExponential(ComplexNumber **A, int nMax, ComplexNumber **res){

ComplexNumber M[3][3];
     for (int i = 0; i< 3; i++ ){
        for (int j = 0; j< 3; j++ ){
            M[i][j] = A[i][j];
        }}
    ComplexNumber zZero;
    ComplexNumber zOne;
    ComplexNumber zTotal;
    ComplexNumber zPrev;
    ComplexNumber zCurrent;
    ComplexNumber zNext;
    ComplexNumber zSum;

for (int i = 0; i< 3; i++ ){
    for (int j = 0; j< 3; j++ ){
        //std::cout << "M00 = " << M[i][j] << "\n";
        zTotal = ComplexNumber (0.0, 0.0); //initialise and to define entries not along diagonals which is just 0
        int factorial = 1; //initialise 0!
        //for a diagonal matrix, exp. of matrix = exp. of each diagonal element
    if (i==j){
        zZero = ComplexNumber(1.0, 0.0); //for n = 0, element always = 1
        zOne = M[i][j]; //for n=1, element always = that matrix entry
        zPrev = M[i][j]; //to initialise the previous one so current n term in series = prev n term*next n term
        zSum = ComplexNumber (0.0, 0.0); //to keep track of sum
    for(int n = 2; n <= nMax; n++)
{
                //std::cout << "Prev = " << zPrev << "\n";
                //to multiply complex numbers
                double a = RealPart(zPrev);
                double b = ImaginaryPart(zPrev);
                double c = RealPart(M[i][j])/double(n);
                double d = ImaginaryPart(M[i][j])/double(n);
                zNext = ComplexNumber(c,d); //next n term in series
                //std::cout << "Next = " << zNext << "\n";
                //multiply prev n term in series with next n term in series
                double Re = (a*c) - (b*d);
                double Im = (c*b) + (a*d);
                zCurrent = ComplexNumber(Re, Im); //then put together again as complex numbers in current n term sum
                //std::cout << "Current = " << zCurrent << "\n";
                zSum = zSum + zCurrent; //add it to the current sum
                zPrev = zCurrent; //prev n term is now the current n sum
                //repeat with multiplying with next n term and adding that to the sum as a fraction (so top and bottom does not blow up)
            }
            zTotal = zZero + zOne + zSum;
            //std::cout << "Total = " << zTotal << "\n";

        res[i][j] = zTotal;}
        else{
            res[i][j] = zTotal;
        }
        //std::cout << "Result = " << res[i][j] << "\n";
    }}
}

/* —— */

/* Derived classes */

//Parent class:
#ifndef SUBMISSION_STUDENT_HPP_
#define SUBMISSION_STUDENT_HPP_

#include <string>
#include <iostream>

class Student {
public:
	Student();
	Student(std::string name, double fines, double fees);

	std::string name;
	double tuition_fees;
	virtual double MoneyOwed() const; //so derived classes can use this and overload it

	void SetLibraryFines(double amount);
	double GetLibraryFines() const;

private:
	double library_fines;
};



#endif /* SUBMISSION_STUDENT_HPP_ */

#include "Student.hpp"
#include <cmath>
#include <iostream>
#include <string>

Student::Student()
{
    name = "unspecified";
    library_fines = 0.0;
    tuition_fees = 0.0;
}
Student::Student(std::string name, double fines, double fees)
{
    this->name = name; //since function parameter needs to be assigned to variable with same name
    //this accesses objects own address and points to class
    tuition_fees = fees;
    library_fines = fines;
}

double Student::MoneyOwed() const
{
    return (library_fines + tuition_fees);
}

//method allows user to set library fines only to non -ve values
void Student::SetLibraryFines(double amount)
{
    library_fines = fabs(amount);
}
//another method to access the private variable libraryfines
double Student::GetLibraryFines() const
{
    return library_fines;
}

//derived class
#ifndef SUBMISSION_GRADUATESTUDENT_HPP_
#define SUBMISSION_GRADUATESTUDENT_HPP_

#include "Student.hpp"

class GraduateStudent : public Student {
private:

public:
	GraduateStudent();
	GraduateStudent(std::string name, double fines, double fees, bool fullTime);
	bool fullTime;
	virtual double MoneyOwed() const;

};



#endif /* SUBMISSION_GRADUATESTUDENT_HPP_ */

#include "GraduateStudent.hpp"
#include <cmath>
#include <iostream>
#include <string>

//use polymorphism to write a method that calculate total money
//owed by a graduate student, who do not pay tuition fees
//will require method for calculating total money owed to be a virtual function
//of parent class
GraduateStudent::GraduateStudent()
{

}
GraduateStudent::GraduateStudent(std::string name, double fines, double fees, bool fullTime):Student(name, fines, fees) //use parent class constructor variable since only fullTime variable is defined in this derived class
{
    fullTime = fullTime;
}
double GraduateStudent::MoneyOwed() const
{
    return Student::GetLibraryFines(); //way to access private members in derived class using functions defined in parent class
}

//derived class of a derived class
#ifndef SUBMISSION_PHDSTUDENT_HPP_
#define SUBMISSION_PHDSTUDENT_HPP_

#include "GraduateStudent.hpp"

class PhdStudent : public GraduateStudent {
public:
	PhdStudent(std::string name, double fines, double fees, bool fullTime);
	virtual double MoneyOwed() const;
};



#endif /* SUBMISSION_PHDSTUDENT_HPP_ */

#include "PhdStudent.hpp"
#include <cmath>
#include <iostream>
#include <string>
PhdStudent::PhdStudent(std::string name, double fines, double fees, bool fullTime):GraduateStudent(name,fines,fees,fullTime) //use the parent class (so derived class of base class) since no variables are defined in this derived class
{

}
double PhdStudent::MoneyOwed() const
{
    return 0.0;
}

/* —— */

#ifndef SUBMISSION_EXERCISE82_HPP_
#define SUBMISSION_EXERCISE82_HPP_

template<typename T>
T CalcAbs(T val) { //function returns same type as type of input parameter val
	// write a single function to calculate absolute value of an int or double
	T result;
	if (val < 0)
    {
        result = -val;
    }
    else
    {
        result = val;
    }
    return result;
}

#ifndef EXCEPTIONDEF
#define EXCEPTIONDEF

#include <string>
class Exception	{
	private:
		std::string mTag, mProblem;
	public:
		Exception(std::string tagString, std::string probString);
		void PrintDebug() const;
};

#endif //EXCEPTIONDEF

#include <iostream>
#include <fstream>
#include "Exception.hpp"

/*The constructors for each of the two classes should take only the probString argument
and should set the tagString member to an appropriate string.
Write a catch block which is able to catch a generic exception
but can also differentiate between these two types of error.*/

//Constructor
Exception::Exception(std::string tagString, std::string probString)
{
    mTag = tagString;
    mProblem = probString;
}

void Exception::PrintDebug () const
{
    std::cerr << "** Error ("<<mTag<<") ** \n";
    std::cerr << "Problem: " << mProblem << "\n\n";
}

#ifndef SUBMISSION_FILENOTOPENEXCEPTION_HPP_
#define SUBMISSION_FILENOTOPENEXCEPTION_HPP_

#include "Exception.hpp"

class FileNotOpenException : public Exception {
	public:
		FileNotOpenException(std::string prob);
};


#endif /* SUBMISSION_FILENOTOPENEXCEPTION_HPP_ */

#include <iostream>
#include <fstream>
#include "FileNotOpenException.hpp"

FileNotOpenException::FileNotOpenException(std::string probString):Exception("FILE NOT OPEN", probString) //kind of like the same as for GraduateStudent and PhDStudent
{

}


#ifndef SUBMISSION_OUTOFRANGEEXCEPTION_HPP_
#define SUBMISSION_OUTOFRANGEEXCEPTION_HPP_

#include "Exception.hpp"

class OutOfRangeException : public Exception {
	public:
	    OutOfRangeException(std::string prob);
};


#endif /* SUBMISSION_OUTOFRANGEEXCEPTION_HPP_ */

#include <iostream>
#include <fstream>
#include "OutOfRangeException.hpp"

OutOfRangeException::OutOfRangeException(std::string probString):Exception("OUT OF RANGE", probString)
{

}



/* ASSIGNMENT 3*/

/* Templates in classes */

//can change a class into a template class that takes any type (int, double, float etc.) as a vector
#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF
#include <cmath>
#include <iostream>
#include <cassert>

template<class T> //<typename T> is also the same, T stands for Type
//actual type of T would be determined by compiler depending on argument passed to the function e.g. Vector<float>
class Vector
{
private:
   T* mData; // data stored in vector //data could be any type
   int mSize; // size of vector //should be an integer since size is always an integer
public:

   Vector(const Vector<T>& otherVector); //from Vector(const Vector& otherVector); since Vector in argument could be any type
   Vector(int size);
   ~Vector(); //destructor
   T GetSize() const; //from int GetSize() const
   T& operator[](int i); //from double& operator[] (int i) //returned type would the same type as the one called in the function, well in this case, the one called in as the class type i.e. float in Vector<float>
   T Read(int i) const; //from double Read(int i) const
   T& operator()(int i); // one-based indexing
   // assignment
   Vector<T>& operator=(const Vector<T>& otherVector); //from Vector& operator=(const Vector& other Vector) since now Vector<T> is the class that should be called & Vector in argument could be any type
   Vector<T> operator+() const; // unary +
   Vector<T> operator-() const; // unary -
   Vector<T> operator+(const Vector<T>& v1) const; // binary +
   Vector<T> operator-(const Vector<T>& v1) const; // binary -
   // scalar multiplication
   Vector<T> operator*(T a) const; //from Vector operator*(double a) since scalar a could be any type
   // p-norm method
   T CalculateNorm(int p=2) const; //from double CalculateNorm(int p=2) const;
   // declare length function as a friend
    //declare all instantiations of template as friends using
    template<typename U> //T is fixed to current class template so need another type for a friend
   friend int length(const Vector<U>& v); //from friend int length(const Vector& v);
};


// Overridden copy constructor
// Allocates memory for new vector, and copies
// entries of other vector into it
template <class T> //always have to call the template first then the class:
Vector<T>::Vector(const Vector<T>& otherVector)
{
   mSize = otherVector.GetSize();
   mData = new T [mSize]; //remember to replace type of mData (which was double) with T
   for (int i=0; i<mSize; i++)
   {
      mData[i] = otherVector.mData[i];
   }
}

// Constructor for vector of a given size
// Allocates memory, and initialises entries
// to zero
template <class T>
Vector<T>::Vector(int size)
{
   assert(size > 0);
   mSize = size;
   mData = new T [mSize]; //remember to replace type of mData (which was double) with T
   for (int i=0; i<mSize; i++)
   {
      mData[i] = 0.0;
   }
}

// Overridden destructor to correctly free memory
template <class T>
Vector<T>::~Vector()
{
   delete[] mData;
}

// Method to get the size of a vector
template <class T>
T Vector<T>::GetSize() const
{
   return mSize;
}

// Overloading square brackets
// Note that this uses `zero-based' indexing,
// and a check on the validity of the index
template <class T>
T& Vector<T>::operator[](int i)
{
   assert(i > -1);
   assert(i < mSize);
   return mData[i];
}

// Read-only variant of []
// Note that this uses `zero-based' indexing,
// and a check on the validity of the index
template <class T>
T Vector<T>::Read(int i) const
{
   assert(i > -1);
   assert(i < mSize);
   return mData[i];
}

// Overloading round brackets
// Note that this uses `one-based' indexing,
// and a check on the validity of the index
template <class T>
T& Vector<T>::operator()(int i)
{
   assert(i > 0);
   assert(i < mSize+1);
   return mData[i-1];
}

// Overloading the assignment operator
template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& otherVector)
{
   assert(mSize == otherVector.mSize);
   for (int i=0; i<mSize; i++)
   {
      mData[i] = otherVector.mData[i];
   }
   return *this;
}

// Overloading the unary + operator
template<class T>
Vector<T> Vector<T>::operator+() const
{
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = mData[i];
   }
   return v;
}
template<class T>
// Overloading the unary - operator
Vector<T> Vector<T>::operator-() const
{
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = -mData[i];
   }
   return v;
}

// Overloading the binary + operator
template<class T>
Vector<T> Vector<T>::operator+(const Vector<T>& v1) const
{
   assert(mSize == v1.mSize);
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = mData[i] + v1.mData[i];
   }
   return v;
}

// Overloading the binary - operator
template<class T>
Vector<T> Vector<T>::operator-(const Vector<T>& v1) const
{
   assert(mSize == v1.mSize);
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = mData[i] - v1.mData[i];
   }
   return v;
}

template<class T>
Vector<T> Vector<T>::operator*(T a) const
{
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = a*mData[i];
   }
   return v;
}

// Method to calculate norm (with default value p=2)
// corresponding to the Euclidean norm
template<class T>
T Vector<T>::CalculateNorm(int p) const
{
   double norm_val, sum = 0.0;
   for (int i=0; i<mSize; i++)
   {
      sum += pow(fabs(mData[i]), p);
   }
   norm_val = pow(sum, 1.0/((double)(p)));
   return norm_val;
}

// MATLAB style friend to get the size of a vector
template<class U> //now need to declare all instantiations as friends of the
int length(const Vector<U>& v)
{
   return v.mSize;
}

#endif
//Code from Chapter10.tex line 19 save as Vector.hpp


//uses Vector.hpp for friend function:
#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF
#include <cmath>
#include <cassert>
#include "Vector.hpp"

template<class M>
class Matrix
{
private:
   M** mData; // entries of matrix //which can be any type
   int mNumRows, mNumCols; // dimensions
public:
   Matrix(const Matrix<M>& otherMatrix);
   Matrix(int numRows, int numCols);
   ~Matrix();
   M GetNumberOfRows() const;
   M GetNumberOfColumns() const;
   M& operator()(int i, int j); //1-based indexing
   //overloaded assignment operator
   Matrix<M>& operator=(const Matrix<M>& otherMatrix);
   Matrix<M> operator+() const; // unary +
   Matrix<M> operator-() const; // unary -
   Matrix<M> operator+(const Matrix<M>& m1) const; // binary +
   Matrix<M> operator-(const Matrix<M>& m1) const; // binary -
   // scalar multiplication
   Matrix<M> operator*(M a) const;
   M CalculateDeterminant() const;
   // declare vector multiplication friendship
   template <class T> //declare the type of class first, T could actually be U as well, not actually declaring the class template here, just changing the variable

   friend Vector<T> operator*(const Matrix<T>& m,
                           const Vector<T>& v);
   template <class T>
   friend Vector<T> operator*(const Vector<T>& v,
                           const Matrix<T>& m);

};


// Overwritten copy constructor
// Allocate memory for new matrix, and copy
// entries into this matrix
template <class M>
Matrix<M>::Matrix(const Matrix<M>& otherMatrix)
{
   mNumRows = otherMatrix.mNumRows;
   mNumCols = otherMatrix.mNumCols;
   mData = new M* [mNumRows];
   for (int i=0; i<mNumRows; i++)
   {
      mData[i] = new M [mNumCols];
   }
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mData[i][j] = otherMatrix.mData[i][j];
      }
   }
}

// Constructor for vector of a given length
// Allocates memory, and initialises entries
// to zero
template<class M>
Matrix<M>::Matrix(int numRows, int numCols)
{
   assert(numRows > 0);
   assert(numCols > 0);
   mNumRows = numRows;
   mNumCols = numCols;
   mData = new M* [mNumRows];
   for (int i=0; i<mNumRows; i++)
   {
      mData[i] = new M [mNumCols];
   }
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mData[i][j] = 0.0;
      }
   }
}

// Overwritten destructor to correctly free memory
template <class M>
Matrix<M>::~Matrix()
{
   for (int i=0; i<mNumRows; i++)
   {
      delete[] mData[i];
   }
   delete[] mData;
}

// Method to get number of rows of matrix
template <class M>
M Matrix<M>::GetNumberOfRows() const
{
   return mNumRows;
}

// Method to get number of columns of matrix
template <class M>
M Matrix<M>::GetNumberOfColumns() const
{
   return mNumCols;
}

// Overloading the round brackets
// Note that this uses `one-based' indexing,
// and a check on the validity of the index
template <class M>
M& Matrix<M>::operator()(int i, int j)
{
   assert(i > 0);
   assert(i < mNumRows+1);
   assert(j > 0);
   assert(j < mNumCols+1);
   return mData[i-1][j-1];
}

// Overloading the assignment operator
template <class M>
Matrix<M>& Matrix<M>::operator=(const Matrix<M>& otherMatrix)
{
   assert(mNumRows = otherMatrix.mNumRows);
   assert(mNumCols = otherMatrix.mNumCols);

   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mData[i][j] = otherMatrix.mData[i][j];
      }
   }
   return *this;
}

// Overloading the unary + operator
template<class M>
Matrix<M> Matrix<M>::operator+() const
{
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = mData[i][j];
      }
   }
   return mat;
}

// Overloading the unary - operator
template <class M>
Matrix<M> Matrix<M>::operator-() const
{
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = -mData[i][j];
      }
   }
   return mat;
}

// Overloading the binary + operator
template <class M>
Matrix<M> Matrix<M>::operator+(const Matrix<M>& m1) const
{
   assert(mNumRows == m1.mNumRows);
   assert(mNumCols == m1.mNumCols);
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = mData[i][j] + m1.mData[i][j];
      }
   }
   return mat;
}

// Overloading the binary - operator
template <class M>
Matrix<M> Matrix<M>::operator-(const Matrix<M>& m1) const
{
   assert(mNumRows == m1.mNumRows);
   assert(mNumCols == m1.mNumCols);
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = mData[i][j] - m1.mData[i][j];
      }
   }
   return mat;
}

// Overloading scalar multiplication
template<class M>
Matrix<M> Matrix<M>::operator*(M a) const
{
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = a*mData[i][j];
      }
   }
   return mat;
}



// Calculate determinant of square matrix recursively
template <class M>
M Matrix<M>::CalculateDeterminant() const
{
   assert(mNumRows == mNumCols);
   M determinant = 0.0;

   if (mNumRows == 1)
   {
      determinant = mData[0][0];
   }
   else
   {
      // More than one entry of matrix
      for (int i_outer=0; i_outer<mNumRows; i_outer++)
      {
         Matrix sub_matrix(mNumRows-1,
                             mNumRows-1);
         for (int i=0; i<mNumRows-1; i++)
         {
            for (int j=0; j<i_outer; j++)
            {
               sub_matrix(i+1,j+1) = mData[i+1][j];
            }
            for (int j=i_outer; j<mNumRows-1; j++)
            {
               sub_matrix(i+1,j+1) = mData[i+1][j+1];
            }
         }
         M sub_matrix_determinant =
                  sub_matrix.CalculateDeterminant();

         determinant += pow(-1.0, i_outer)*
                  mData[0][i_outer]*sub_matrix_determinant;
      }
   }
   return determinant;
}

    // Overloading matrix multiplied by a vector
    template <class T>
    Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v)
    {
       int original_vector_size = v.GetSize();
       assert(m.GetNumberOfColumns() == original_vector_size);
       int new_vector_length = m.GetNumberOfRows();
       Vector<T> new_vector(new_vector_length);

       for (int i=0; i<new_vector_length; i++)
       {
          for (int j=0; j<original_vector_size; j++)
          {
             new_vector[i] += m.mData[i][j]*v.Read(j);
          }
       }

       return new_vector;
    }
   // Overloading vector multiplied by a matrix
   template <class T>
    Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m)
    {
       int original_vector_size = v.GetSize();
       assert(m.GetNumberOfRows() == original_vector_size);
       int new_vector_length = m.GetNumberOfColumns();
       Vector<T> new_vector(new_vector_length);

       for (int i=0; i<new_vector_length; i++)
       {
          for (int j=0; j<original_vector_size; j++)
          {
             new_vector[i] += v.Read(j)*m.mData[j][i];
          }
       }

       return new_vector;
}



#endif
//Code from Appendix.tex line 608 save as Matrix.hpp


/* ASSIGNMENT 4*/

//using std::vector to represent 1D array of Matrix elements and represent them in row or column major

#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF
#include <cmath>
#include <cassert>
#include <vector>
#include "Vector.hpp"

template<class M>
class Matrix
{
private:
   int mNumRows, mNumCols; // dimensions
   std::vector<M> mData;
   bool isRowMajor;
public:
   Matrix(const Matrix<M>& otherMatrix);
   Matrix(int numRows, int numCols, bool rowMajor);
   ~Matrix();
   M GetNumberOfRows() const;
   M GetNumberOfColumns() const;
   M& operator()(int i, int j); //1-based indexing
   //overloaded assignment operator
   Matrix<M>& operator=(const Matrix<M>& otherMatrix);
   Matrix<M> operator+() const; // unary +
   Matrix<M> operator-() const; // unary -
   Matrix<M> operator+(const Matrix<M>& m1) const; // binary +
   Matrix<M> operator-(const Matrix<M>& m1) const; // binary -
   // scalar multiplication
   Matrix<M> operator*(M a) const;
   M CalculateDeterminant() const;
   // declare vector multiplication friendship
   template <class T>
   friend Vector<T> operator*(const Matrix<T>& m,
                           const Vector<T>& v);
   template <class T>
   friend Vector<T> operator*(const Vector<T>& v,
                           const Matrix<T>& m);
    //new public functions in Matrix.hpp
    std::vector<M> getData() const;
    bool getOrder() const;


};


// Overwritten copy constructor
// Allocate memory for new matrix, and copy
// entries into this matrix
template <class M>
Matrix<M>::Matrix(const Matrix<M>& otherMatrix)
{
   mNumRows = otherMatrix.mNumRows;
   mNumCols = otherMatrix.mNumCols;
   isRowMajor = otherMatrix.isRowMajor;
   mData.resize(mNumRows*mNumCols);
   for (int i = 0; i < (mNumRows*mNumCols); i++)
   {
       mData[i] = otherMatrix.mData[i];
   }

}

// Constructor for vector of a given length
// Allocates memory, and initialises entries
// to zero
template<class M>
Matrix<M>::Matrix(int numRows, int numCols, bool rowMajor)
{
   assert(numRows > 0);
   assert(numCols > 0);
   mNumRows = numRows;
   mNumCols = numCols;
   isRowMajor = rowMajor;
   mData.resize(mNumRows*mNumCols);
   for (int i = 0; i < mNumRows; i++)
   {
       for (int j = 0; j < mNumCols; j++)
       {
          mData[i*j] = 0.0;
      }
   }
}

// Overwritten destructor to correctly free memory
template <class M>
Matrix<M>::~Matrix()
{

}


// Method to get number of rows of matrix
template <class M>
M Matrix<M>::GetNumberOfRows() const
{
   return mNumRows;
}

// Method to get number of columns of matrix
template <class M>
M Matrix<M>::GetNumberOfColumns() const
{
   return mNumCols;
}

// Overloading the round brackets
// Note that this uses `one-based' indexing,
// and a check on the validity of the index
template <class M>
M& Matrix<M>::operator()(int i, int j)
{
   assert(i > 0);
   assert(i < mNumRows+1);
   assert(j > 0);
   assert(j < mNumCols+1);
   //change this part
   if (isRowMajor == 1)
   {
       return mData[((i-1)*mNumCols) + (j-1)];
   }
   else
   {
       return mData[(i-1) + ((j-1)*mNumRows)];
   }
}

// Overloading the assignment operator
template <class M>
Matrix<M>& Matrix<M>::operator=(const Matrix<M>& otherMatrix)
{
   assert(mNumRows = otherMatrix.mNumRows);
   assert(mNumCols = otherMatrix.mNumCols);
    mData.resize(mNumRows*mNumCols);

   for (int i = 0; i < (mNumRows*mNumCols); i++)
   {
       mData[i] = otherMatrix.mData[i];
   }

   return *this;
}

// Overloading the unary + operator
template<class M>
Matrix<M> Matrix<M>::operator+() const
{
   Matrix mat(mNumRows, mNumCols, isRowMajor);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
        if (isRowMajor == 1)
            {mat(i+1,j+1) = mData[(i*mNumCols) + j];}
        else
            {mat(i+1,j+1) = mData[i + (j*mNumRows)];}
      }}
   return mat;

}

// Overloading the unary - operator
template <class M>
Matrix<M> Matrix<M>::operator-() const
{
   Matrix mat(mNumRows, mNumCols, isRowMajor);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
        if (isRowMajor == 1)
            {mat(i+1,j+1) = -mData[(i*mNumCols) + j];}
        else
            {mat(i+1,j+1) = -mData[i + (j*mNumRows)];}
      }}
   return mat;

}

// Overloading the binary + operator
template <class M>
Matrix<M> Matrix<M>::operator+(const Matrix<M>& m1) const
{
   assert(mNumRows == m1.mNumRows);
   assert(mNumCols == m1.mNumCols);
   assert(isRowMajor == m1.isRowMajor);
   Matrix mat(mNumRows, mNumCols, isRowMajor);
  for (int i=0; i<mNumRows; i++)
    {
  for (int j=0; j<mNumCols; j++)
    {
    if (isRowMajor == 1)
        {mat(i+1,j+1) = mData[(i*mNumCols) + j] + m1.mData[(i*mNumCols) + j];}
    else
        {mat(i+1,j+1) = mData[i + (j*mNumRows)] + m1.mData[i + (j*mNumRows)];}
    }}
   return mat;
}

// Overloading the binary - operator
template <class M>
Matrix<M> Matrix<M>::operator-(const Matrix<M>& m1) const
{
   assert(mNumRows == m1.mNumRows);
   assert(mNumCols == m1.mNumCols);
   assert(isRowMajor == m1.isRowMajor);
   Matrix mat(mNumRows, mNumCols, isRowMajor);
  for (int i=0; i<mNumRows; i++)
    {
  for (int j=0; j<mNumCols; j++)
    {
    if (isRowMajor == 1)
        {mat(i+1,j+1) = mData[(i*mNumCols) + j] - m1.mData[(i*mNumCols) + j];}
    else
        {mat(i+1,j+1) = mData[i + (j*mNumRows)] - m1.mData[i + (j*mNumRows)];}
    }}
   return mat;
}

// Overloading scalar multiplication
template<class M>
Matrix<M> Matrix<M>::operator*(M a) const
{
   Matrix mat(mNumRows, mNumCols, isRowMajor);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
        if (isRowMajor == 1)
            {mat(i+1,j+1) = a*mData[(i*mNumCols) + j];}
        else
            {mat(i+1,j+1) = a*mData[i + (j*mNumRows)];}
      }}
   return mat;

}



// Calculate determinant of square matrix recursively
template <class M>
M Matrix<M>::CalculateDeterminant() const
{
   assert(mNumRows == mNumCols);
   M determinant = 0.0;

   if (mNumRows == 1)
   {
      determinant = mData[0];
   }
   else
   {
      // More than one entry of matrix
      for (int i_outer=0; i_outer<mNumRows; i_outer++)
      {
         Matrix sub_matrix(mNumRows-1,
                             mNumRows-1);
         for (int i=0; i<mNumRows-1; i++)
         {
            for (int j=0; j<i_outer; j++)
            {
                if (isRowMajor == 1)
                {sub_matrix(i+1,j+1) = mData[((i+1)*mNumCols) + j];}
                else
                {sub_matrix(i+1,j+1) = mData[(i+1) + (j*mNumRows)];}
            }
            for (int j=i_outer; j<mNumRows-1; j++)
            { if (isRowMajor == 1)
                {sub_matrix(i+1,j+1) = mData[((i+1)*mNumCols) + (j+1)];}
               else
                {sub_matrix(i+1,j+1) = mData[(i+1) + ((j+1)*mNumRows)];}
            }
         }
         M sub_matrix_determinant =
                  sub_matrix.CalculateDeterminant();

         determinant += pow(-1.0, i_outer)*
                  mData[0][i_outer]*sub_matrix_determinant;
      }
   }
   return determinant;
}

    // Overloading matrix multiplied by a vector
    template <class T>
    Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v)
    {
       int original_vector_size = v.GetSize();
       assert(m.GetNumberOfColumns() == original_vector_size);
       int new_vector_length = m.GetNumberOfRows();
       bool rowMajor = m.getOrder();
       int numCols = m.GetNumberOfColumns();
       int numRows = m.GetNumberOfRows();
       Vector<T> new_vector(new_vector_length);

       for (int i=0; i<new_vector_length; i++)
       {
          for (int j=0; j<original_vector_size; j++)
          {
              if (rowMajor == 1)
            {new_vector[i] += m.mData[(i*numCols) + j]*v.Read(j);}
                else
            {new_vector[i] += m.mData[i + (j*numRows)]*v.Read(j);}
          }
       }

       return new_vector;
    }
   // Overloading vector multiplied by a matrix
   template <class T>
    Vector<T> operator*(const Vector<T>& v, const Matrix<T>& m)
    {
       int original_vector_size = v.GetSize();
       assert(m.GetNumberOfRows() == original_vector_size);
       int new_vector_length = m.GetNumberOfColumns();
        bool rowMajor = m.getOrder();
       int numCols = m.GetNumberOfColumns();
       int numRows = m.GetNumberOfRows();
       Vector<T> new_vector(new_vector_length);

       for (int i=0; i<new_vector_length; i++)
       {
          for (int j=0; j<original_vector_size; j++)
          {
              if (rowMajor == 1)
            {new_vector[i] += v.Read(j)*m.mData[(j*numCols) + i];}
                else
            {new_vector[i] += v.Read(j)*m.mData[j + (i*numRows)];}
          }
       }

       return new_vector;
}

template <class M>
std::vector<M> Matrix<M>::getData() const
{
    return mData;
}

template <class M>
bool Matrix<M>::getOrder() const
{
    return isRowMajor;
}


#endif
//Code from Appendix.tex line 608 save as Matrix.hpp


