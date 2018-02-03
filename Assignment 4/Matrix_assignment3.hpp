#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF
#include <cmath>
#include <cassert>
#include "Vector.hpp"

template<class M>
class Matrix
{
private:
   M** mData; // entries of matrix
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
   template <class T>
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
