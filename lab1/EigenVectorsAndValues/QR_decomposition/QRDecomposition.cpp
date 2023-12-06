#include "QRDecomposition.h"

QRDecomposition::QRDecomposition(Matrix &A) :
    matrix(A),
    rows(A.rows),
    cols(A.cols)
{
    this->Q = Matrix(rows, cols);
    this->R = Matrix(rows, cols);
}

std::vector<std::string> QRDecomposition::find(double epsilon, int maxIterations)
{
    if (rows != cols)
    {
        throw std::runtime_error("Expected square matrix, but rectangular found! Therefore inversion of given matrix is impossible.");
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (std::abs(matrix.data[i][j]) < epsilon) {
                matrix.data[i][j] = epsilon;
            }
        }
    }

    for (; iterations < maxIterations; ++iterations)
    {
        if (eigenValuesReady(matrix, epsilon)) break;
        qr_decomposition();
        matrix = R * Q;
    }

//    std::ofstream out("../FilesWithResults/outfile.txt");
//    out << "Матрица:" << std::endl << matrix << std::endl;
//    out.close();

    if (iterations == maxIterations - 1) throw std::runtime_error("Max iterations");

    std::vector<std::string> eigenValues(rows);

    int i = 0;
    for (; i < rows - 1; ++i)
    {
        if (std::abs(matrix.data[i + 1][i]) < epsilon) {
            eigenValues[i] = std::to_string(matrix.data[i][i]);
            continue;
        }

        double b = -matrix.data[i][i] - matrix.data[i + 1][i + 1];
        double c = matrix.data[i][i] * matrix.data[i + 1][i + 1] - matrix.data[i][i + 1] * matrix.data[i + 1][i];
        double d = b * b - 4 * c;

        if (d >= 0)
        {
            d = std::sqrt(d);
            eigenValues[i] = std::to_string( -b / 2.0 + d / 2.0);
            eigenValues[i + 1] = std::to_string( -b / 2.0 - d / 2.0);
            i += 1;
            continue;
        }

        d = std::sqrt(-d);
        eigenValues[i] = std::to_string(-b / 2.0) + "+" + std::to_string(d / 2.0) + "i";
        eigenValues[i + 1] = std::to_string(-b / 2.0) + "-" + std::to_string(d / 2.0) + "i";
        i += 1;
    }

    if (i == rows - 1) eigenValues[i] = std::to_string(matrix.data[i][i]);

    return eigenValues;
}

void QRDecomposition::qr_decomposition()
{
    if (rows != cols)
    {
        throw std::runtime_error("Expected square matrix, but rectangular found! Therefore inversion of given matrix is impossible.");
    }

    R = matrix;
    Q = Matrix(rows, cols);
    for (int i = 0; i < rows; ++i) Q.data[i][i] = 1.0;

    for (int i = 0; i < cols; ++i)
    {
        Matrix V(rows, 1);
        double element = R.data[i][i];

        int sign = 1;
        if (element < 0.0) sign = -1;

        double sum = element * element;
        for (int j = i + 1; j < rows; ++j)
        {
            V.data[j][0] = R.data[j][i];
            sum += R.data[j][i] * R.data[j][i];
        }

        V.data[i][0] = element + sign * std::sqrt(sum);

        Matrix E(rows, cols);
        for (int j = 0; j < rows; ++j) E.data[j][j] = 1.0;

        double div = (~V * V).data[0][0];
        if (div < 1e-10)
        {
            div = 1e-10;
        }

        Matrix H = E - V * ~V * (2.0 / div);

        R = H * R;
        Q = Q * ~H;
    }
}

bool QRDecomposition::eigenValuesReady(const Matrix &A, double epsilon)
{
    static int passed_blocks = 0;
    static int passed_blocks_save = 0;

    for (int i = 0; i < A.cols - 1; ++i)
    {
        double sum = 0.0;
        double element = A.data[i + 1][i] * A.data[i + 1][i];
        for (int j = i + 2; j < A.rows; ++j)
        {
            sum += A.data[j][i] * A.data[j][i];
        }
        if (std::sqrt(sum + element) < epsilon) continue;

        if (std::sqrt(sum) < element && element >= epsilon)
        {
            double b = -A.data[i][i] - A.data[i + 1][i + 1];
            double c = A.data[i][i] * A.data[i + 1][i + 1] - A.data[i][i + 1] * A.data[i + 1][i];
            double d = b * b - 4 * c;

            static std::pair<double, double> previous_lambda(0,0);

            if (d >= 0)
            {
                d = std::sqrt(d);
                std::pair<double, double> delta( -b / 2 + d / 2 - previous_lambda.first, 0 );
                previous_lambda = { -b / 2 + d / 2, 0 };

                if (std::abs(delta.first) < epsilon)
                {
                    previous_lambda = { 0, 0 };
                    passed_blocks_save += 1;
                    i += 1;
                    continue;
                }
                passed_blocks = passed_blocks_save;
                return false;
            }

            if (passed_blocks > 0)
            {
                passed_blocks -= 1;
                i += 1;
                continue;
            }

            d = std::sqrt(-d);
            std::pair<double, double> delta( -b / 2.0 - previous_lambda.first, d / 2.0 - previous_lambda.second );
            previous_lambda = { -b / 2.0, d / 2.0 };

            if (std::sqrt(delta.first * delta.first + delta.second * delta.second) < epsilon)
            {
                previous_lambda = { 0, 0 };
                passed_blocks_save += 1;
                i += 1;
                continue;
            }

            passed_blocks = passed_blocks_save;
            return false;
        }

        passed_blocks = passed_blocks_save;
        return false;
    }

    return true;
}

Matrix QRDecomposition::get_R() const
{
    return R;
}
Matrix QRDecomposition::get_Q() const
{
    return Q;
}

long long QRDecomposition::getCountIterations() const
{
    return iterations;
}