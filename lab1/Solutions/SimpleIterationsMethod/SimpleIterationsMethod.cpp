#include "SimpleIterationsMethod.h"

SimpleIterations::SimpleIterations(Matrix &A) :
    AbstractSolution(A)
{

}

std::vector<double> SimpleIterations::solution(std::vector<double> vector, double epsilon, int maxIterations)
{
    if (!isCorrectMatrix())
    {
        throw std::runtime_error("Incorrect matrix!");
    }

    std::vector<double> X(rows);

    //Начальное приближение
    for (int i = 0; i < rows; ++i)
    {
        X[i] = vector[i] / matrix.data[i][i];
    }

    for (; iterations < maxIterations; ++iterations)
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

bool SimpleIterations::isCorrectMatrix() const
{
    for (int i = 0; i < rows; ++i)
    {
        double sum = 0;
        for (int j = i + 1; j < cols; ++j) sum += fabs(matrix.data[i][j]);
        if (fabs(matrix.data[i][i]) < fabs(sum)) return false;
    }

    return true;
}

long long SimpleIterations::getCountIterations() const
{
    return iterations;
}