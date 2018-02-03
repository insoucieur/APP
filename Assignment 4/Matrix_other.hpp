/**
 * 2D Matrix with non-contiguous 2D storage internally.
 */
template <typename TValue>
class Matrix
{
public:
   /**
    * the value type.
    */
   typedef TValue ValueType;

   /**
    * the value type.
    */
   typedef Matrix<ValueType> Self;

   /**
    * the pointer type.
    */
   typedef ValueType* PointerType;

   /**
    * the storage type.
    */
   typedef ValueType** StorageType;

   /**
    * the size type.
    */
   typedef int SizeType;

   /**
    * contructor.
    */
   Matrix( SizeType nrow, SizeType ncol )
      : m_nrow( nrow )
      , m_ncol( ncol )
      , m_data( Create( nrow, ncol ) )
   {
   }

   /**
    * destructor.
    */
   virtual ~Matrix()
   {
      for ( int row = 0; row < Rows(); ++row )
      {
         delete[] m_data[ row ];
      }

      delete[] m_data;
   }

   /**
    * the i-th row;
    * this allows you to write matrix[y][x] to access an element.
    */
   PointerType operator[]( int i )
   {
      return m_data[ i ];
   }

   /**
    * alias for Get(y,x) (const);
    * this allows you to write int x = matrix(y,x) to access an element.
    */
   ValueType operator()( SizeType y, SizeType x ) const
   {
      return Get( y, x );
   }

   /**
    * alias for Get(y,x) (non-const);
    * this allows you to write matrix(y,x) = x; to access an element.
    */
   ValueType& operator()( SizeType y, SizeType x )
   {
      return Get( y, x );
   }

   /**
    * the number of rows.
    */
   SizeType Rows() const
   {
      return m_nrow;
   }

   /**
    * the number of columns.
    */
   SizeType Columns() const
   {
      return m_ncol;
   }

   /**
    * the value of the given location (const).
    */
   ValueType Get( SizeType y, SizeType x ) const
   {
      return m_data[ y ][ x ];
   }

   /**
    * the value of the given location (non-const).
    */
   ValueType& Get( SizeType y, SizeType x )
   {
      return m_data[ y ][ x ];
   }

   /**
    * set the value of the given location.
    */
   void Set( SizeType y, SizeType x, const ValueType& value )
   {
      m_data[ y ][ x ] = value;
   }

   /**
    * the internal representation.
    */
   friend StorageType GetImpl( Self& matrix )
   {
      return matrix.m_data;
   }

protected:
   /**
    * convenience function to create the internal representation inside the
    * member initialization list (MIL).
    */
   StorageType Create( SizeType nrow, SizeType ncol )
   {
      StorageType m = new PointerType[ nrow ];

      for ( int row = 0; row < nrow; ++row )
      {
         m[ row ] = new ValueType[ ncol ];
      }

      return m;
   }

private:
   /**
    * the number of columns.
    */
   SizeType m_ncol;

   /**
    * the number of rows.
    */
   SizeType m_nrow;

   /**
    * the internal representation.
    */
   StorageType m_data;

private:
   /**
    * prevent asignement.
    */
   Self& operator= ( const Self& );

   /**
    * prevent copy-asignement.
    */
   Matrix( const Self& );
};
