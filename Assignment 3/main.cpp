#include "Vector.hpp"
#include "Matrix.hpp"

int main ()
{
    Vector<float> v(2); //must explicitly pass the template type in class templates
    v(1) = 3.14;

    Matrix<int> m(2,2);
    m(2,2) = 13;

    return 0;
}
