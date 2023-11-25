#include <iostream>
#include <tuple>
#include "Solutions/LU_decomposition/LUDecomposition.h"
#include "Solutions/SweepMethod/SweepMethod.h"
#include "Solutions/SimpleIterationsMethod/SimpleIterationsMethod.h"
#include "Solutions/SeidelMethod/SeidelMethod.h"
#include "EigenVectorsAndValues/RotationMethod/RotationMethod.h"
#include "EigenVectorsAndValues/QR_decomposition/QRDecomposition.h"

//Вывод вектора
std::ostream& operator<<(std::ostream &out, const std::vector<double> &vector)
{
    for (auto element : vector) out << element << std::endl;
    return out;
}

//Считывание вектора
std::istream& operator>>(std::istream& in, std::vector<double> &vector)
{
    for (double &element : vector) in >> element;
    return in;
}

//Получение данных с файла для нахождения СВ и СЗ
std::tuple<Matrix, double, int> getResourcesForEigenVectors(const std::string &fileName)
{
    std::ifstream inputFile(fileName);
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

//Получение данных с файла для решения СЛАУ
std::tuple<Matrix, std::vector<double>> getResources(const std::string fileName)
{
    std::ifstream inputFile(fileName);
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

//Вывод исходных данных для СВ и СЗ
void outputOfInitialValues(std::ostream &out, const Matrix &A)
{
    out << "Исходная матрица:" << std::endl;
    out << A << std::endl;
}

//Вывод исходных данных для СЛАУ
void outputOfInitialValues(std::ostream &out, const Matrix &A, const std::vector<double> &vector)
{
    out << "Исходная матрица:" << std::endl;
    out << A << std::endl;
    out << "Вектор столбец:" << std::endl;
    out << vector << std::endl;
}

//Нахождение СВ и СЗ
std::tuple<Matrix, std::vector<double>> findEigenVectorsAndValues(AbstractFinder &finder, const double epsilon, const int maxIterations)
{
    return finder.find(epsilon, maxIterations);
}

//Решение СЛАУ
std::vector<double> solveSystem(AbstractSolution &solver, const std::vector<double> &vector, const double epsilon = 0.1, const int maxIterations = 1000)
{
    return solver.solution(vector, epsilon, maxIterations);
}

//Вывод найденных СВ и СЗ
void outputEigenVectorsAndValues(std::ostream &out, const std::tuple<Matrix, std::vector<double>> &answer)
{
    auto EigenVectors = std::get<0>(answer);
    auto EigenValues = std::get<1>(answer);

    out << "Собственные векторы:" << std::endl;
    for (int i = 0; i < EigenVectors.getRows(); ++i)
    {
        out << "Вектор #" << i + 1 << ": ";
        EigenVectors.printCol(out, i);
        out << std::endl;
    }

    out << std::endl << "Собственные значения:" << std::endl;
    for (int i = 0; i < EigenValues.size(); ++i) out << "Значение " << i + 1 << " = " << EigenValues[i] << std::endl;
}

//Вывод решения СЛАУ
void outputOfTheSystemSolution(std::ostream &out, const std::vector<double> &vector)
{
    out << "Решение СЛАУ:" << std::endl;
    for (int i = 0; i < vector.size(); ++i) out << "x" << i + 1 << " = " << vector[i] << std::endl;
}

//Проверка найденных СВ и СЗ
void check(std::ostream &out, const Matrix &A, const std::tuple<Matrix, std::vector<double>> &answer)
{
    auto EigenVectors = std::get<0>(answer);
    auto EigenValues = std::get<1>(answer);

    for (int i = 0; i < EigenVectors.getCols(); ++i)
    {
        std::vector<double> EigenVector = EigenVectors.getCol(i);

        std::vector<double> multiplyMatrixAndEigenVector = A * EigenVector;

        out << "Для собственного вектора (";
        for (int j = 0; j < EigenVector.size(); ++j) out << EigenVector[j] << ((j != EigenVector.size() - 1) ? ", " : ")^T ");
        out << "и собственного значения " << EigenValues[i] << std::endl;

        out << "A * EigenVector" << std::endl << "(";
        for (int j = 0; j < EigenVector.size(); ++j) out << multiplyMatrixAndEigenVector[i] << ((j != EigenVector.size() - 1) ? ", " : ")^T");

        out << std::endl;

        std::vector<double> multiplyEigenValueAndEigenVector(EigenVector.size());
        for (int j = 0; j < EigenVector.size(); ++j) multiplyEigenValueAndEigenVector[j] = EigenValues[i] * EigenVector[j];

        out << "EigenValue * EigenVector" << std::endl << "(";
        for (int j = 0; j < EigenVector.size(); ++j) out << multiplyEigenValueAndEigenVector[i] << ((j != EigenVector.size() - 1) ? ", " : ")^T");

        out << std::endl << std::endl;
    }
}

void task1(const std::string &inputFileName, const std::string &outputFileName)
{
    auto resources = getResources(inputFileName);
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

void task2(const std::string &inputFileName, const std::string &outputFileName)
{
    auto resources = getResources(inputFileName);
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

void task3(const std::string &inputFileName, const std::string &outputFileName)
{
    auto resources = getResources(inputFileName);
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

    double eps = 0.001;
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

void task4(const std::string &inputFileName, const std::string &outputFileName)
{
    auto resources = getResourcesForEigenVectors(inputFileName);
    Matrix A = std::get<0>(resources);
    double epsilon = std::get<1>(resources);
    int maxIterations = std::get<2>(resources);

    std::ofstream outputFile(outputFileName);

    outputOfInitialValues(outputFile, A);

    RotationMethod finder(A);

    try
    {
        outputFile << "Собственных векторы и значения, найденные методом вращений с точностью epsilon = " << epsilon << std::endl;
        auto answer = findEigenVectorsAndValues(finder, epsilon, maxIterations);
        outputEigenVectorsAndValues(outputFile, answer);
        outputFile << std::endl << "Количество итераций " << finder.getCountIterations() << std::endl;
        outputFile << std::endl << "Проверка:" << std::endl;
        check(outputFile, A, answer);
        std::cout << "Done! Check file outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
    }

    outputFile.close();
}

void task5(const std::string &inputFileName, const std::string &outputFileName)
{
    auto resources = getResourcesForEigenVectors(inputFileName);
    Matrix A = std::get<0>(resources);
    double epsilon = std::get<1>(resources);
    int maxIterations = std::get<2>(resources);

    std::ofstream outputFile(outputFileName);

    outputOfInitialValues(outputFile, A);

    QRDecomposition finder(A);

    try
    {
        outputFile << "Собственных векторы и значения, найденные с помощью QR разложения с точностью epsilon = " << epsilon << std::endl;
        auto answer = finder.find(epsilon, maxIterations);
        for (int i = 0; i < answer.size(); ++i) outputFile << "Значение #" << i + 1 << " = " << answer[i] << std::endl;
        outputFile << std::endl << "Количество итераций " << finder.getCountIterations() << std::endl;
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

    const std::string inputFileName = "../FilesWithResults/input.txt";
    const std::string outputFileName = "../FilesWithResults/outfile.txt";

//    task1(inputFileName, outputFileName);
//    task2(inputFileName, outputFileName);
//    task3(inputFileName, outputFileName);
//    task4(inputFileName, outputFileName);
    task5(inputFileName, outputFileName);

    return 0;
}

