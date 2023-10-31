#include "Matrix.h"

Matrix::Matrix() :
        rows(0), cols(0),
        data(0, std::vector<double>(0))
{

}

Matrix::Matrix(int n, int m) :
        rows(n),
        cols(m),
        data(n, std::vector<double>(m))
{

}

void Matrix::inputMatrix()
{
    std::cout << "Enter matrix elements:" << std::endl;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            std::cin >> data[i][j];
        }
    }
}

void Matrix::inputMatrixFromFile(std::ifstream &filename)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            filename >> data[i][j];
        }
    }
}

void Matrix::transform_l(Matrix &L, int i, int j)
{
    for (int k = 0; k < i; ++k)
    {
        std::swap(L.data[i][k], L.data[j][k]);
    }
}

void Matrix::swap_row(Matrix & matrix, int i, int j)
{
    std::swap(matrix.data[i], matrix.data[j]);
}

int Matrix::find_pivot_row(const Matrix & matrix, int col)
{
    int pivot = col;
    for (int i = col + 1; i < matrix.rows; i++)
    {
        if (fabs(matrix.data[i][col]) > fabs(matrix.data[pivot][col]))
        {
            pivot = i;
        }
    }

    return pivot;
}

std::ostream& operator<<(std::ostream &out, const Matrix& matrix)
{
    for (int i = 0; i < matrix.rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            out << matrix.data[i][j] << " ";
        }
        out << std::endl;
    }
    return out;
}

Matrix operator+(const Matrix &A, const Matrix &B)
{
    if (A.rows != B.rows || A.cols != B.cols)
    {
        throw;
    }

    Matrix C(A.rows, A.cols);

    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            C.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }

    return C;
}

Matrix operator-(const Matrix &A, const Matrix &B)
{
    if (A.rows != B.rows || A.cols != B.cols)
    {
        throw;
    }

    Matrix C(A.rows, A.cols);

    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            C.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }

    return C;
}

std::vector<double> operator*(const std::vector<double> &vector, const Matrix &B)
{
    return B * vector;
}

std::vector<double> operator*(const Matrix &A, const std::vector<double> &vector)
{
    if (A.cols != vector.size())
    {
        throw;
    }

    std::vector<double> C(A.rows);

    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            C[i] += A.data[i][j] * vector[j];
        }
    }

    return C;
}

Matrix operator*(int x, const Matrix &A)
{
    return A * x;
}

Matrix operator*(const Matrix &A, int x)
{
    Matrix C(A.rows, A.cols);

    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            C.data[i][j] = A.data[i][j] * x;
        }
    }

    return C;
}

Matrix operator*(const Matrix &A, const Matrix &B)
{
    if (A.cols != B.rows)
    {
        throw;
    }

    Matrix C(A.rows, B.cols);

    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < B.cols; ++j)
        {
            for (int k = 0; k < A.cols; ++k)
            {
                C.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }

    return C;
}

void Matrix::operator*() const
{
    std::cout << *this;
}

Matrix operator~(const Matrix &A)
{
    Matrix C(A.cols, A.rows);

    for (int i = 0; i < C.rows; ++i)
    {
        for(int j = 0; j < C.cols; ++j)
        {
            C.data[i][j] = A.data[j][i];
        }
    }

    return C;
}
