#include "Matrix.h"

Matrix::Matrix(int n, int m) :
        rows(n),
        cols(m),
        output_file("../output.txt"),
        debug_file("../debug.txt"),
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

std::vector<double> Matrix::solve(std::vector<double> b) const
{
    if (get_determinant() == 0)
    {
        throw std::runtime_error("Given matrix has zero or infinity solutions!");
    }

    Matrix L(this->rows, this->cols);
    Matrix U(this->rows, this->cols);
    Matrix P(this->rows, this->cols);

    int not_used = 0;

    this->lu_decomposition(L, U, P, not_used);

    b = ~P * b;

    std::vector<double> Z(L.rows);

    for (int i = 0; i < L.rows; ++i)
    {
        double sum = 0;
        for (int j = 0; j < i; ++j)
        {
            sum += L.data[i][j] * Z[j];
        }

        Z[i] = b[i] - sum;
    }

    std::vector<double> X(U.rows);

    for (int i = U.rows -1; i >= 0; --i)
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

void Matrix::lu_decomposition(Matrix &L, Matrix &U, Matrix &P, int &count) const
{
    if (this->rows != this->cols)
    {
        return;
    }

    U = *this;

    for (int i = 0; i < U.rows; i++)
    {
        L.data[i][i] = 1.0;
        P.data[i][i] = 1.0;
    }

    std::ofstream file_for_debug(this->debug_file);

    file_for_debug << "Исходные значения" << std::endl << std::endl;
    file_for_debug << "Матрица L:" << std::endl << L << std::endl;
    file_for_debug << "Матрица U:" << std::endl << U << std::endl;
    file_for_debug << "Матрица P:" << std::endl << P;

    int iteration = 1;

    for (int i = 0; i < U.rows - 1; i++)
    {
        int pivot = find_pivot_row(U, i);
        if (pivot != i)
        {
            swap_row(U, i, pivot);
            swap_row(P, i, pivot);
            transform_l(L, i, pivot);
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

        file_for_debug << "---------------------------" << std::endl;
        file_for_debug << "Итеация №" << iteration++ << std::endl << std::endl;
        file_for_debug << "Матрица L:" << std::endl << L << std::endl;
        file_for_debug << "Матрица U:" << std::endl << U << std::endl;
        file_for_debug << "Матрица P:" << std::endl << ~P;
    }

    P = ~P;
}

double Matrix::get_determinant() const
{
    if (this->rows != this->cols)
    {
        throw std::runtime_error("Cannot find determinant of non square matrix!");
    }

    Matrix L(this->rows, this->cols);
    Matrix U(this->rows, this->cols);
    Matrix P(this->rows, this->cols);

    int count_changes = 0;

    this->lu_decomposition(L, U, P, count_changes);

    double determinant(1.0);

    for (int i = 0; i < U.rows; ++i)
    {
        determinant *= U.data[i][i];
    }

    return (count_changes & 1) ? -determinant : determinant;
}

Matrix Matrix::get_inverse() const
{
    if (this->rows != this->cols)
    {
        throw std::runtime_error("Expected square matrix, but rectangular found! Therefore inversion of given matrix is impossible.");
    }

    if (this->get_determinant() == 0)
    {
        throw std::runtime_error("Determinant of matrix equals zero thus its not inversable!");
    }

    Matrix identity_matrix(this->rows, this->cols);

    for (int i = 0; i < this->rows; ++i)
    {
        identity_matrix.data[i][i] = 1.0;
    }

    Matrix L(this->rows, this->cols);
    Matrix U(this->rows, this->cols);
    Matrix P(this->rows, this->cols);
    int not_used = 0;

    lu_decomposition(L, U, P, not_used);

    std::vector<double> B(this->rows);

    Matrix inverse_matrix(this->rows, this->cols);

    for (int i = 0; i < this->rows; ++i)
    {
        for (int j = 0; j < this->cols; ++j)
        {
            B[j] = identity_matrix.data[j][i];
        }

        std::vector<double> partial_solution = solve(B);

        for (int k = 0; k < this->rows; ++k)
        {
            inverse_matrix.data[k][i] = partial_solution[k];
        }
    }

    return inverse_matrix;
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
