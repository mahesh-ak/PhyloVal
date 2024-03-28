#ifndef MATRIX_H
#define MATRIX_H



#include<iostream>

#include"include_blas.hpp"



template<typename E>
class Matrix {
private:

    int height, width;
    unsigned long allocatedSize;
    E* data;



public:

    typedef CBLAS_TRANSPOSE Transpose;
    static const Transpose TRANSPOSE = CblasTrans;
    static const Transpose NO_TRANSPOSE = CblasNoTrans;

    Matrix()
        : height(1), width(1),
          allocatedSize(100),
          data(new E[allocatedSize])
    {}

    Matrix(int height, int width)
        : height(height), width(width),
          allocatedSize((unsigned long)height * width),
          data(new E[allocatedSize])
    {}

    ~Matrix() {
        delete[] data;
    }



private:

    // Not implemented, dont copy matrices.
    Matrix(const Matrix &copy);
    Matrix &operator=(const Matrix &copy);

    unsigned Idx(int i, int j) const { return i * width + j; }



public:

    E  operator()(int i, int j) const { return data[Idx(i, j)]; }
    E &operator()(int i, int j)       { return data[Idx(i, j)]; }



    void resize(int height, int width)
    {
        this->height = height;
        this->width = width;
        unsigned long newAllocated = height * width;
        if(newAllocated > allocatedSize)
        {
            allocatedSize = newAllocated;
            delete[] data;
            data = new E[allocatedSize];
        }
    }

    int GetHeight() const { return height; }
    int GetWidth()  const { return width;  }



    static void Mult(const Matrix &in1, Transpose transpose1,
                     const Matrix &in2, Transpose transpose2,
                     Matrix &out)
    {
        cblas_dgemm(CblasRowMajor,
                    transpose1,
                    transpose2,
                    out.height, out.width, transpose1 == NO_TRANSPOSE ? in1.width : in1.height,
                    1.0,
                    in1.data, in1.width,
                    in2.data, in2.width,
                    0.0,
                    out.data, out.width);
    }

    static void Mult(const Matrix &in1,
                     const Matrix &in2,
                     Matrix &out)
    { Mult(in1, NO_TRANSPOSE, in2, NO_TRANSPOSE, out); }
};



template<typename E>
inline
std::ostream &operator<<(std::ostream &out, const Matrix<E> &m)
{
    out << "[" << m.GetHeight() << "," << m.GetWidth() << "]";
    out << "(";
    std::string sep1 = "";
    for(int i = 0; i < m.GetHeight(); ++i)
    {
        out << sep1 << "(";
        std::string sep2 = "";
        for(int j = 0; j < m.GetWidth(); ++j)
        {
            out << sep2 << m(i, j);
            sep2 = ",";
        }
        out << ")";
        sep1 = ",";
    }
    out << ")";

    return out;
}



#endif
