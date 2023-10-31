#include "LU_decomposition.h"

LU_decomposition::LU_decomposition(Matrix &A)
{
    _A = A;
    Matrix L(A.rows, A.cols);
    Matrix U(A.rows, A.cols);
    Matrix P(A.rows, A.cols);
    int count = 0;
    lu_decomposition(A, L, U, P, count);
    _L = L;
    _U = U;
    _P = P;
    _count_permutations = count;
}

void LU_decomposition::lu_decomposition(Matrix &A, Matrix &L, Matrix &U, Matrix &P, int &count)
{
    if (A.rows != A.cols)
    {
        return;
    }

    U = A;

    for (int i = 0; i < U.rows; i++)
    {
        L.data[i][i] = 1.0;
        P.data[i][i] = 1.0;
    }

    for (int i = 0; i < U.rows - 1; i++)
    {
        int pivot = Matrix::find_pivot_row(U, i);
        if (pivot != i)
        {
            Matrix::swap_row(U, i, pivot);
            Matrix::swap_row(P, i, pivot);
            Matrix::transform_l(L, i, pivot);
            count++;
        }

        for (int j = i + 1; j < U.rows; j++)
        {
            double factor = U.data[j][i] / U.data[i][i];
            L.data[j][i] = factor;

            for (int k = i; k < U.rows; k++)
            {
                U.data[j][k] -= U.data[i][k] * factor;
            }
        }
    }

    P = ~P;
}

std::vector<double> LU_decomposition::solve(std::vector<double> b)
{
    if (get_determinant() == 0)
    {
        throw std::runtime_error("Given matrix has zero or infinity solutions!");
    }

    b = ~_P * b;

    std::vector<double> Z(_L.rows);

    for (int i = 0; i < _L.rows; ++i)
    {
        double sum = 0;
        for (int j = 0; j < i; ++j)
        {
            sum += _L.data[i][j] * Z[j];
        }

        Z[i] = b[i] - sum;
    }

    std::vector<double> X(_U.rows);

    for (int i = _U.rows -1; i >= 0; --i)
    {
        double sum = 0;
        for (int j = _U.cols; j >= 1; --j)
        {
            sum += _U.data[i][j] * X[j];
        }

        X[i] = (Z[i] - sum) / _U.data[i][i];
    }

    return X;
}

double LU_decomposition::get_determinant()
{
    if (_A.rows != _A.cols)
    {
        throw std::runtime_error("Cannot find determinant of non square matrix!");
    }

    double determinant(1.0);

    for (int i = 0; i < _U.rows; ++i)
    {
        determinant *= _U.data[i][i];
    }

    return (_count_permutations & 1) ? -determinant : determinant;
}

Matrix LU_decomposition::get_inverse()
{
    if (_A.rows != _A.cols)
    {
        throw std::runtime_error("Expected square matrix, but rectangular found! Therefore inversion of given matrix is impossible.");
    }

    if (this->get_determinant() == 0)
    {
        throw std::runtime_error("Determinant of matrix equals zero thus its not inversable!");
    }

    Matrix identity_matrix(_A.rows, _A.cols);

    for (int i = 0; i < _A.rows; ++i)
    {
        identity_matrix.data[i][i] = 1.0;
    }

    std::vector<double> B(_A.rows);

    Matrix inverse_matrix(_A.rows, _A.cols);

    for (int i = 0; i < _A.rows; ++i)
    {
        for (int j = 0; j < _A.cols; ++j)
        {
            B[j] = identity_matrix.data[j][i];
        }

        std::vector<double> partial_solution = solve(B);

        for (int k = 0; k < _A.rows; ++k)
        {
            inverse_matrix.data[k][i] = partial_solution[k];
        }
    }

    return inverse_matrix;
}

Matrix LU_decomposition::get_L() const
{
    return _L;
}

Matrix LU_decomposition::get_U() const
{
    return _U;
}

Matrix LU_decomposition::get_P() const
{
    return _P;
}

Matrix LU_decomposition::get_PLU() const
{
    Matrix result = _P * _L * _U;
    return result;
}