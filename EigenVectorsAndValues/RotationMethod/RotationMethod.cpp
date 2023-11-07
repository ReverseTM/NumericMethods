#include "RotationMethod.h"

RotationMethod::RotationMethod(Matrix &A) :
    AbstractFinder(A)
{

}

std::tuple<Matrix, std::vector<double>> RotationMethod::find(double epsilon, int maxIterations)
{
    if (!isCorrectMatrix())
    {
        throw std::runtime_error("Incorrect matrix");
    }

    Matrix A = matrix;

    Matrix V(rows, cols);
    for (int i = 0; i < rows; ++i) V.data[i][i] = 1.0;

    for (; iterations < maxIterations; ++iterations)
    {
        double maxValue = 0.0;
        double sum = 0.0;
        int p = 0, q = 0;

        for (int i = 0; i < rows; ++i)
        {
            for (int j = i + 1; j < cols; ++j)
            {
                sum += A.data[i][j] * A.data[i][j];
                if (fabs(A.data[i][j]) > maxValue)
                {
                    maxValue = fabs(A.data[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        double stop_iteration = sqrt(sum);

        if (stop_iteration < epsilon) break;

        double angle = 0.5 * atan(2 * A.data[p][q] / (A.data[p][p] - A.data[q][q]));
//                (A.data[p][p] == A.data[q][q])
//                ? (M_PI / 4)
                //: 0.5 * atan(2 * A.data[p][q] / (A.data[p][p] - A.data[q][q]));

        Matrix U(rows, cols);
        for (int i = 0; i < rows; ++i) U.data[i][i] = 1.0;

        U.data[p][p] = cos(angle);
        U.data[q][q] = cos(angle);
        U.data[p][q] = -sin(angle);
        U.data[q][p] = sin(angle);

        A = (~U * A) * U;
        V = V * U;
    }

    std::vector<double> eigenValues(rows);
    for (int i = 0; i < rows; ++i) eigenValues[i] = A.data[i][i];

    return std::make_tuple(V, eigenValues);
}

bool RotationMethod::isCorrectMatrix() const
{
    return (matrix == ~matrix);
}

long long RotationMethod::getCountIterations() const
{
    return iterations;
}