#include "Matrix.hpp"

#include <iostream>
#include <cstdlib>
#include <cmath>



int main(int argn, char** argv)
{
    Matrix<double> m1(2, 1);
    m1(0, 0) = 1;
    m1(1, 0) = 2;

    Matrix<double> m2(1, 3);
    m2(0, 0) = 2;
    m2(0, 1) = 3;
    m2(0, 2) = 5;

    Matrix<double> out(2, 3);

    Matrix<double>::Mult(m1, m2, out);

    for(unsigned i = 0; i < 2; ++i)
    {
        for(unsigned j = 0; j < 3; ++j)
        {
            if(fabs(out(i, j) - m1(i, 0) * m2(0, j)) > 0.0000000001)
            {
                std::cout << "Entry at (" << i << ", " << j << ") is wrong." << std::endl;
                std::exit(-1);
            }
        }
    }



    Matrix<double> out2(3, 3);
    Matrix<double>::Mult(m2, Matrix<double>::TRANSPOSE,
                         m2, Matrix<double>::NO_TRANSPOSE,
                         out2);

    for(unsigned i = 0; i < 3; ++i)
    {
        for(unsigned j = 0; j < 3; ++j)
        {
            if(fabs(out2(i, j) - m2(0, i) * m2(0, j)) > 0.0000000001)
            {
                std::cout << "Entry at (" << i << ", " << j << ") is wrong." << std::endl;
                std::exit(-1);
            }
        }
    }

    return 0;
}
