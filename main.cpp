#include <iostream>
#include <tuple>
#include "Solutions/LU_decomposition/LUDecomposition.h"
#include "Solutions/SweepMethod/SweepMethod.h"
#include "Solutions/SimpleIterationsMethod/SimpleIterationsMethod.h"
#include "Solutions/SeidelMethod/SeidelMethod.h"
#include "EigenVectorsAndValues/RotationMethod/RotationMethod.h"

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

std::tuple<Matrix, double, int> getResourcesForEigenVectors()
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

    double epsilon;
    inputFile >> epsilon;

    int maxIterations;
    inputFile >> maxIterations;

    inputFile.close();

    return std::make_tuple(A, epsilon, maxIterations);
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

void outputOfInitialValues(std::ostream &out, Matrix &A)
{
    out << "Исходная матрица:" << std::endl;
    out << A << std::endl;
}

void outputOfInitialValues(std::ostream &out, Matrix &A, std::vector<double> &vector)
{
    out << "Исходная матрица:" << std::endl;
    out << A << std::endl;
    out << "Вектор столбец:" << std::endl;
    out << vector << std::endl;
}

std::vector<double> solveSystem(AbstractSolution &solver, std::vector<double> &vector, double epsilon = 0.1, int maxIterations = 1000)
{
    return solver.solution(vector, epsilon, maxIterations);
}

void outputEigenVectorsAndValues(std::ostream &out, std::tuple<Matrix, std::vector<double>> &answer)
{
    auto EigenVectors = std::get<0>(answer);
    auto EigenValues = std::get<1>(answer);

    out << "Собственные векторы:" << std::endl;
    for (int i = 0; i < EigenVectors.getRows(); ++i)
    {
        out << "Вектор #" << i + 1 << ": ";
        EigenVectors.printRow(out, i);
        out << std::endl;
    }

    out << std::endl << "Собственные значения:" << std::endl;
    for (int i = 0; i < EigenValues.size(); ++i) out << "Значение " << i + 1 << " = " << EigenValues[i] << std::endl;
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

    try
    {
        std::vector<double> answer = solveSystem(solver, vector);
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
        outputFile.close();
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
        outputFile.close();
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

    try
    {
        std::vector<double> answer = solveSystem(solver, vector);
        outputOfTheSystemSolution(outputFile, answer);
        std::cout << "Done! Check file outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
    }

    outputFile.close();
}

void task3()
{
    auto resources = getResources();
    Matrix A = std::get<0>(resources);
    std::vector vector = std::get<1>(resources);

    std::ofstream outputFile(outputFileName);

    outputOfInitialValues(outputFile, A, vector);

    std::vector<AbstractSolution*> solvers
    {
        new SimpleIterationsMethod(A),
        new SeidelMethod(A)
    };

    std::vector<std::string> methodName = {"Метод простых итераций", "Метод Зейделя"};

    double eps = 0.000001;
    int maxIterations = 1000;

    try
    {
        for (int i = 0; i < solvers.size(); ++i)
        {
            outputFile << "Решение используя " + methodName[i] << " с точностью epsilon = " << eps << std::endl;
            std::vector<double> answer = solveSystem(*(solvers[i]), vector, eps, maxIterations);
            outputOfTheSystemSolution(outputFile, answer);
            outputFile << "Количество итераций " << ((SeidelMethod*)solvers[i])->getCountIterations() << std::endl << std::endl;
        }
        std::cout << "Успешно! Проверьте файл outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
        for (auto solver : solvers) delete solver;
    }

    outputFile.close();
    for (auto solver : solvers) delete solver;
}

void task4()
{
    auto resources = getResourcesForEigenVectors();
    Matrix A = std::get<0>(resources);
    double epsilon = std::get<1>(resources);
    int maxIterations = std::get<2>(resources);

    std::ofstream outputFile(outputFileName);

    outputOfInitialValues(outputFile, A);

    RotationMethod solver(A);

    try
    {
        outputFile << "Значения собственных векторов и значений, найденных методом вращений с точностью epsilon = " << epsilon << std::endl;
        auto answer = solver.find(epsilon, maxIterations);
        outputEigenVectorsAndValues(outputFile, answer);
        outputFile << std::endl << "Количество итераций " << solver.getCountIterations() << std::endl;
        std::cout << "Done! Check file outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
    }

    outputFile.close();
}

int main() {

    task4();
    return 0;
}

