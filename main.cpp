#include <iostream>
#include "Matrix/Matrix.h"

void task1()
{
    std::string filename = "../LU_decomposition/input.txt";
    std::ifstream input_file(filename);

    int n;
    input_file >> n;

    Matrix A(n, n);

    A.inputMatrixFromFile(input_file);

    std::vector<double> B(n);

    for (int i = 0; i < n; ++i)
    {
        input_file >> B[i];
    }

    input_file.close();

    std::ofstream out_file("../outfile.txt");

    out_file << "Исходная матрица матрица:" << std::endl << A << std::endl;
    out_file << "Вектор столбец:" << std::endl;
    for (auto item : B)
    {
        out_file << item << std::endl;
    }

    out_file << std::endl;

    out_file << "LU разложение:" << std::endl;

    Matrix L(n, n);
    Matrix U(n, n);
    Matrix P(n, n);
    int not_used = 0;
    A.lu_decomposition(A, L, U, P, not_used);

    out_file << "Матрица L:" << std::endl << L << std::endl;
    out_file << "Матрица U:" << std::endl << U << std::endl;
    out_file << "Матрица P:" << std::endl << P << std::endl;
    out_file << "Произведение P * L * U:" << std::endl << P * L * U << std::endl;


    out_file << "Решение СЛАУ:" << std::endl;
    try
    {
        std::vector<double> result = A.solve(B);
        for (int i = 0; i < B.size(); ++i)
        {
            out_file << "x" << i + 1 << " = " << result[i] << std::endl;
        }

        out_file << std::endl;
    }
    catch (std::exception &ex)
    {
        out_file << ex.what();
    }


    out_file << "Определитель матрицы:" << std::endl;
    try
    {
        out_file << A.get_determinant() << std::endl;

        out_file << std::endl;
    }
    catch (std::exception &ex)
    {
        out_file << ex.what();
    }

    out_file << "Обратная матрица:" << std::endl;

    try
    {
        auto inverse = A.get_inverse();
        out_file << inverse << std::endl;

        out_file << "Произведение исходной матрицы на обратную A * A^-1 = E:" << std::endl;
        out_file << A * inverse << std::endl;

    }
    catch (std::exception &ex)
    {
        out_file << ex.what();
    }

    out_file.close();

    std::cout << "Done! Check file outfile.txt" << std::endl;
}

int main() {


    return 0;
}

