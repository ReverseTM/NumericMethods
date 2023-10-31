#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(Matrix &A)
{
    this->matrix = A;
    this->rows = A.rows;
    this->cols = A.cols;
    this->L = Matrix(rows, cols);
    this->U = Matrix(rows, cols);
    this->P = Matrix(rows, cols);

    lu_decomposition();
}

std::vector<double> LUDecomposition::solution(std::vector<double> vector)
{
    if (get_determinant() == 0)
    {
        throw std::runtime_error("Given matrix has zero or infinity solutions!");
    }

    vector = ~P * vector;

    std::vector<double> Z(rows);

    for (int i = 0; i < rows; ++i)
    {
        double sum = 0;
        for (int j = 0; j < i; ++j)
        {
            sum += L.data[i][j] * Z[j];
        }

        Z[i] = vector[i] - sum;
    }

    std::vector<double> X(rows);

    for (int i = rows -1; i >= 0; --i)
    {
        double sum = 0;
        for (int j = U.cols; j >= 1; --j)
        {
            sum += U.data[i][j] * X[j];
        }

        X[i] = (Z[i] - sum) / U.data[i][i];
    }

    return X;
}

void LUDecomposition::lu_decomposition()
{
    if (rows != cols)
    {
        return;
    }

    U = this->matrix;

    for (int i = 0; i < rows; i++)
    {
        L.data[i][i] = 1.0;
        P.data[i][i] = 1.0;
    }

    for (int i = 0; i < rows - 1; i++)
    {
        int pivot = Matrix::find_pivot_row(U, i);
        if (pivot != i)
        {
            Matrix::swap_row(U, i, pivot);
            Matrix::swap_row(P, i, pivot);
            Matrix::transform_l(L, i, pivot);
            count_permutations++;
        }

        for (int j = i + 1; j < rows; j++)
        {
            double factor = U.data[j][i] / U.data[i][i];
            L.data[j][i] = factor;

            for (int k = i; k < rows; k++)
            {
                U.data[j][k] -= U.data[i][k] * factor;
            }
        }
    }

    P = ~P;
}

double LUDecomposition::get_determinant()
{
    if (rows != cols)
    {
        throw std::runtime_error("Cannot find determinant of non square matrix!");
    }

    double determinant(1.0);

    for (int i = 0; i < U.rows; ++i)
    {
        determinant *= U.data[i][i];
    }

    return (count_permutations & 1) ? -determinant : determinant;
}

Matrix LUDecomposition::get_inverse()
{
    if (rows != cols)
    {
        throw std::runtime_error("Expected square matrix, but rectangular found! Therefore inversion of given matrix is impossible.");
    }

    if (this->get_determinant() == 0)
    {
        throw std::runtime_error("Determinant of matrix equals zero!");
    }

    Matrix identity_matrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        identity_matrix.data[i][i] = 1.0;
    }

    std::vector<double> B(rows);

    Matrix inverse_matrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            B[j] = identity_matrix.data[j][i];
        }

        std::vector<double> partial_solution = solution(B);

        for (int k = 0; k < rows; ++k)
        {
            inverse_matrix.data[k][i] = partial_solution[k];
        }
    }

    return inverse_matrix;
}

Matrix LUDecomposition::get_L() const
{
    return L;
}

Matrix LUDecomposition::get_U() const
{
    return U;
}

Matrix LUDecomposition::get_P() const
{
    return P;
}

Matrix LUDecomposition::get_PLU() const
{
    Matrix result = P * L * U;
    return result;
}