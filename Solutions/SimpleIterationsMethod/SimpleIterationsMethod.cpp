#include "SimpleIterationsMethod.h"

SimpleIterationsMethod::SimpleIterationsMethod(Matrix &A, double eps)
{
    matrix = A;
    rows = A.rows;
    cols = A.cols;
    epsilon = eps;
}

std::vector<double> SimpleIterationsMethod::solution(std::vector<double> vector)
{
    if (!isCorrectMatrix())
    {
        throw std::runtime_error("Incorrect matrix!");
    }

    std::vector<double> X(rows);

    int max_iterations = 1000;

    //Начальное приближение
    for (int i = 0; i < rows; ++i)
    {
        X[i] = vector[i] / matrix.data[i][i];
    }

    for (; iterations < max_iterations; ++iterations)
    {
        //Следующие приближения
        std::vector<double> Xn(rows);

        for (int i = 0; i < rows; ++i)
        {
            double sum_term = 0;
            for (int j = 0; j < cols; ++j)
            {
                if (i == j) continue;

                sum_term += matrix.data[i][j] * X[j];
            }

            Xn[i] = (vector[i] - sum_term) / matrix.data[i][i];
        }

        bool flag = true;
        for (int i = 0; i < rows - 1; ++i)
        {
            if (fabs(Xn[i] - X[i]) > epsilon)
            {
                flag = false;
                break;
            }
        }

        X = Xn;

        if (flag) return X;
    }

    throw std::runtime_error("Max iterations");
}

bool SimpleIterationsMethod::isCorrectMatrix() const
{
    for (int i = 0; i < rows; ++i)
    {
        double sum = 0;
        for (int j = i + 1; j < cols; ++j) sum += fabs(matrix.data[i][j]);
        if (fabs(matrix.data[i][i]) < fabs(sum)) return false;
    }

    return true;
}

long long SimpleIterationsMethod::getCountIterations() const
{
    return iterations;
}