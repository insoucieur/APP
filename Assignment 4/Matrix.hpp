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
