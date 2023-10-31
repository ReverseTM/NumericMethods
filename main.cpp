#include <iostream>
#include "Solutions/LU_decomposition/LUDecomposition.h"
#include "Solutions/Sweep_method/SweepMethod.h"

std::string inputFileName = "../FilesWithResults/input.txt";
std::string outputFileName = "../FilesWithResults/outfile.txt";

void task1()
{
    std::ifstream input_file(inputFileName);

    if (!input_file.is_open())
    {
        std::cout << "File not open!" << std::endl;
        return;
    }

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

    std::ofstream out_file(outputFileName);

    out_file << "Исходная матрица:" << std::endl << A << std::endl;
    out_file << "Вектор столбец:" << std::endl;
    for (auto item : B)
    {
        out_file << item << std::endl;
    }

    out_file << std::endl;

    out_file << "LU разложение:" << std::endl;

    LUDecomposition solver(A);

    out_file << "Матрица L:" << std::endl << solver.get_L() << std::endl;
    out_file << "Матрица U:" << std::endl << solver.get_U() << std::endl;
    out_file << "Матрица P:" << std::endl << solver.get_P() << std::endl;
    out_file << "Произведение P * L * U:" << std::endl << solver.get_PLU() << std::endl;


    out_file << "Решение СЛАУ:" << std::endl;
    try
    {
        auto result = solver.solution(B);
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
        out_file << solver.get_determinant() << std::endl;

        out_file << std::endl;
    }
    catch (std::exception &ex)
    {
        out_file << ex.what();
    }

    out_file << "Обратная матрица:" << std::endl;

    try
    {
        auto inverse = solver.get_inverse();
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

void task2()
{
    std::ifstream input_file(inputFileName);

    if (!input_file.is_open())
    {
        std::cout << "File not open!" << std::endl;
        return;
    }

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

    SweepMethod solver(A);

    std::ofstream out_file(outputFileName);

    out_file << "Исходная матрица:" << std::endl << A << std::endl;
    out_file << "Вектор столбец:" << std::endl;
    for (auto item : B)
    {
        out_file << item << std::endl;
    }

    out_file << std::endl;

    out_file << "Решение СЛАУ:" << std::endl;
    try
    {
        auto result = solver.solution(B);
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

    out_file.close();

    std::cout << "Done! Check file outfile.txt" << std::endl;
}

int main() {

    task2();
    return 0;
}

