#include "SweepMethod.h"

SweepMethod::SweepMethod(Matrix &A) :
    AbstractSolution(A),
    P(A.rows),
    Q(A.rows)
{

}

std::vector<double> SweepMethod::solution(std::vector<double> vector, double not_used1, int not_used2)
{
    if (!is_correct_matrix())
    {
        throw std::runtime_error("Matrix is not correct!");
    }

    // Для первой строки
    P[0] = (-matrix.data[0][1]) / matrix.data[0][0];
    Q[0] = vector[0] / matrix.data[0][0];

    // Прямой ход
    for (int i = 1; i < rows - 1; ++i)
    {
        P[i] = (-matrix.data[i][i + 1]) / (matrix.data[i][i] + matrix.data[i][i - 1] * P[i - 1]);
        Q[i] = (vector[i] - matrix.data[i][i - 1] * Q[i - 1]) / (matrix.data[i][i] + matrix.data[i][i - 1] * P[i - 1]);
    }

    // Для последней строки
    P[rows - 1] = 0;
    Q[rows - 1] =
            (vector[rows - 1] - matrix.data[rows - 1][rows - 2] * Q[rows - 2])
            /
            (matrix.data[rows - 1][rows - 1] + matrix.data[rows - 1][rows - 2] * P[rows - 2]);

    std::vector<double> result(rows);

    result[rows - 1] = Q[rows - 1];
    for (int i = rows - 1; i > 0; --i)
    {
        result[i - 1] = P[i - 1] * result[i] + Q[i - 1];
    }

    return result;
}

bool SweepMethod::is_correct_matrix() const
{
    if (
         fabs(matrix.data[0][0]) < fabs(matrix.data[0][1])
         ||
         fabs(matrix.data[rows - 1][rows - 1]) < fabs(matrix.data[rows - 1][rows - 2])
       )
    {
        return false;
    }

    if (matrix.data[0][0] == 0)
    {
        return false;
    }

    for (int i = 1; i < rows - 1; ++i)
    {
        if (matrix.data[i][i] == 0)
        {
            return false;
        }
        if (fabs(matrix.data[i][i]) < fabs(matrix.data[i][i - 1]) + fabs(matrix.data[i][i + 1]))
        {
            return false;
        }
    }

    if (matrix.data[rows - 1][rows - 1] == 0)
    {
        return false;
    }

    return true;
}

std::vector<double> SweepMethod::get_P() const
{
    return P;
}
std::vector<double> SweepMethod::get_Q() const
{
    return Q;
}