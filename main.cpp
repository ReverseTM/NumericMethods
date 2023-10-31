#include <iostream>
#include <tuple>
#include "Solutions/LU_decomposition/LUDecomposition.h"
#include "Solutions/Sweep_method/SweepMethod.h"

std::string inputFileName = "../FilesWithResults/input.txt";
std::string outputFileName = "../FilesWithResults/outfile.txt";

std::ostream& operator<<(std::ostream &out, const std::vector<double> &vector)
{
    for (auto element : vector) out << element << std::endl;
    return out;
}

std::istream& operator>>(std::istream& in, std::vector<double> &vector)
{
    for (double &element : vector) in >> element;
    return in;
}

std::tuple<Matrix, std::vector<double>> getResources()
{
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open())
    {
        std::cout << "File not open!" << std::endl;
        exit(0);
    }

    int n;
    inputFile >> n;

    Matrix A(n, n);
    A.inputMatrixFromFile(inputFile);

    std::vector<double> B(n);
    inputFile >> B;

    inputFile.close();

    return std::make_tuple(A, B);
}

void outputOfInitialValues(std::ostream &out, Matrix &A, std::vector<double> &vector)
{
    out << "Исходная матрица:" << std::endl;
    out << A << std::endl;
    out << "Вектор столбец:" << std::endl;
    out << vector << std::endl;
}

std::vector<double> solveSystem(AbstractSolution &solver, std::vector<double> &vector)
{
    return solver.solution(vector);
}

void outputOfTheSystemSolution(std::ostream &out, std::vector<double> &vector)
{
    out << "Решение СЛАУ:" << std::endl;
    for (int i = 0; i < vector.size(); ++i) out << "x" << i + 1 << " = " << vector[i] << std::endl;
}

void task1()
{
    auto resources = getResources();
    Matrix A = std::get<0>(resources);
    std::vector vector = std::get<1>(resources);

    std::ofstream outputFile(outputFileName);

    outputOfInitialValues(outputFile, A, vector);

    LUDecomposition solver(A);

    std::vector<double> answer;
    try
    {
        answer = solveSystem(solver, vector);
        outputOfTheSystemSolution(outputFile, answer);
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
        return;
    }

    outputFile << "LU разложение:" << std::endl;
    outputFile << "Матрица L:" << std::endl << solver.get_L() << std::endl;
    outputFile << "Матрица U:" << std::endl << solver.get_U() << std::endl;
    outputFile << "Матрица P:" << std::endl << solver.get_P() << std::endl;
    outputFile << "Произведение P * L * U:" << std::endl << solver.get_PLU() << std::endl;

    outputFile << std::endl << "Определитель матрицы:" << std::endl;
    try
    {
        outputFile << solver.get_determinant() << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
    }

    outputFile << std::endl <<  "Обратная матрица:" << std::endl;
    try
    {
        auto inverse = solver.get_inverse();
        outputFile << inverse << std::endl;

        outputFile << "Произведение исходной матрицы на обратную A * A^-1 = E:" << std::endl;
        outputFile << A * inverse << std::endl;

    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
    }

    outputFile.close();

    std::cout << "Done! Check file outfile.txt" << std::endl;
}

void task2()
{
    auto resources = getResources();
    Matrix A = std::get<0>(resources);
    std::vector vector = std::get<1>(resources);

    std::ofstream outputFile(outputFileName);

    outputOfInitialValues(outputFile, A, vector);

    SweepMethod solver(A);

    std::vector<double> answer;
    try
    {
        answer = solveSystem(solver, vector);
        outputOfTheSystemSolution(outputFile, answer);
        std::cout << "Done! Check file outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
    }

    outputFile.close();
}

int main() {

    task1();
    return 0;
}

