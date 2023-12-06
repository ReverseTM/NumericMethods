#include <iostream>
#include <tuple>
#include "lab1/Solutions/LU_decomposition/LUDecomposition.h"
#include "lab1/Solutions/SweepMethod/SweepMethod.h"
#include "lab1/Solutions/SimpleIterationsMethod/SimpleIterationsMethod.h"
#include "lab1/Solutions/SeidelMethod/SeidelMethod.h"
#include "lab1/EigenVectorsAndValues/RotationMethod/RotationMethod.h"
#include "lab1/EigenVectorsAndValues/QR_decomposition/QRDecomposition.h"
#include "lab2/SimpleIterationsMethod/SimpleIterations.h"
#include "lab2/NewtonMethod/NewtonMethod.h"
#include "lab2/SimpleIterationsMethod/SISystem.h"
#include "lab2/NewtonMethod/NMSystem.h"
#include "lab3/Polynomials/LagrangePolynomial.h"
#include "lab3/Polynomials/NewtonPolynomial.h"
#include "lab3/Differentiation/Differentiation.h"
#include "lab3/Integrals/AbstractIntegrater.h"
#include "lab3/Integrals/Rectangle/RectangleMethod.h"
#include "lab3/Integrals/Trapezoid/TrapezoidMethod.h"
#include "lab3/Integrals/Simpson/SimpsonMethod.h"


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
        new SimpleIterations(A),
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

void task6(const std::string &inputFileName, const std::string &outputFileName)
{
    double x0 = 0.0;
    double epsilon = 0.001;
    int maxIterations = 10000;

    auto myFunction { [](double x) {return x * exp(x) + x * x - 1;} };
    auto derivative { [](double x) {return exp(x) + x * exp(x) + 2 * x;} };

    std::ofstream outputFile(outputFileName);

    try
    {
        double x;
        int iterations;

        auto result1 = SimpleIterationsMethod::solve(
                [](double x) {return (x - 0.1345 * ( x * exp(x) + x * x - 1 ));},
                x0,
                epsilon,
                maxIterations
                );

        x = std::get<0>(result1);
        iterations = std::get<1>(result1);

        outputFile << "С помощью метода простых итераций найден корень = " << x << " за " << iterations << " итераций с точностью " << epsilon << std::endl;

        outputFile << "Подставив найденный x в исходное уравнение: F(" << x << ") = " << myFunction(x) << std::endl;

        auto result2 = NewtonMethod::solve(
                myFunction,
                derivative,
                x0,
                epsilon,
                maxIterations
                );

        x = std::get<0>(result2);
        iterations = std::get<1>(result2);

        outputFile << "С помощью метода Ньютона найден корень = " << x << " за " << iterations << " итераций с точностью " << epsilon << std::endl;

        outputFile << "Подставив найденный x в исходное уравнение: F(" << x << ") = " << myFunction(x) << std::endl;

        std::cout << "Done! Check file outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
    }

    outputFile.close();

}

void task7(const std::string &inputFileName, const std::string &outputFileName)
{
    double x0 = 0.0;
    double epsilon = 0.00001;
    int maxIterations = 10000;

    std::vector<double (*)(std::vector<double>)> functionsForSimpleIterations =
    {
            [](std::vector<double> x){ return cos(x[1]) / 2; },
            [](std::vector<double> x) { return exp(x[0]) / 2; }
    };

    auto func1 = [](std::vector<double> x) { return 2 * x[0] - cos(x[1]); };
    auto func2 = [](std::vector<double> x) { return 2 * x[1] - exp(x[0]); };

    std::vector<double (*)(std::vector<double>)> functionsForNewtonMethod =
    {
            func1,
            func2
    };

    std::vector<double (*)(std::vector<double>)> derivativesForNewtonMethod =
    {
            [](std::vector<double> x) { return 2.0; }, //  df1/dx1
            [](std::vector<double> x) { return sin(x[1]); }, //  df1/dx2
            [](std::vector<double> x) { return -exp(x[0]); }, //  df2/dx1
            [](std::vector<double> x) { return 2.0; } //  df2/dx2
    };

    std::ofstream outputFile(outputFileName);

    try
    {
        std::vector<double> x;
        int iterations;

        auto result1 = SISystem::solve(
                functionsForSimpleIterations,
                x0,
                epsilon,
                maxIterations
        );

        x = std::get<0>(result1);
        iterations = std::get<1>(result1);

        outputFile << "С помощью метода простых итераций с точностью " << epsilon << " найдены корни:" << std::endl;
        for (int i = 0; i < x.size(); ++i) outputFile << "x" << i + 1 << " = " << x[i] << std::endl;
        outputFile << "совершенно " << iterations << " итераций." << std::endl;

        outputFile << "Подставив найденные x в исходные уравнения:" << std::endl;
        outputFile << "F1(" << x[0] << ", " << x[1] << ") = " << func1(x) << std::endl;
        outputFile << "F2(" << x[0] << ", " << x[1] << ") = " << func2(x) << std::endl;

        auto result2 = NMSystem::solve(
                functionsForNewtonMethod,
                derivativesForNewtonMethod,
                x0,
                epsilon,
                maxIterations
        );

        x = std::get<0>(result2);
        iterations = std::get<1>(result2);
        outputFile << "---------------------------------------------------------------------" << std::endl;
        outputFile << "С помощью метода Ньютона с точностью " << epsilon << " найдены корни:" << std::endl;
        for (int i = 0; i < x.size(); ++i) outputFile << "x" << i + 1 << " = " << x[i] << std::endl;
        outputFile << "совершенно " << iterations << " итераций." << std::endl;

        outputFile << "Подставив найденные x в исходные уравнения:" << std::endl;
        outputFile << "F1(" << x[0] << ", " << x[1] << ") = " << func1(x) << std::endl;
        outputFile << "F2(" << x[0] << ", " << x[1] << ") = " << func2(x) << std::endl;

        std::cout << "Done! Check file outfile.txt" << std::endl;
    }
    catch (std::exception &ex)
    {
        outputFile << ex.what();
        outputFile.close();
    }

    outputFile.close();
}

void task8(const std::string &inputFileName, const std::string &outputFileName)
{
    std::ofstream out(outputFileName);

    std::vector<double> x2 = {0.1, 0.5, 0.9, 1.3, 1.5};
    std::vector<double> x1 = {0.1, 0.5, 1.1, 1.3, 1.5};
    double x0 = 0.8;

    auto function = [](double x) { return log(x) + x; };

    std::string polynomialString;
    double result;

    auto answer1 = NewtonPolynomial::interpolation(
            function,
            x1,
            x0);

    polynomialString = std::get<0>(answer1);
    result = std::get<1>(answer1);

    out << polynomialString << std::endl;
    out << "P" + std::to_string(x1.size() - 1) + "(" << x0 << ") = " << result << std::endl;
    out << "f(" << x0 << ") = " << function(x0) << std::endl;
    out << "Погрешность " << result - function(x0) << std::endl;

    out << "-------------------------------------------" << std::endl;

    auto answer2 = LagrangePolynomial::interpolation(
            function,
            x2,
            x0);


    polynomialString = std::get<0>(answer2);
    result = std::get<1>(answer2);

    out << polynomialString << std::endl;
    out << "L" + std::to_string(x2.size() - 1) + "(" << x0 << ") = " << result << std::endl;
    out << "f(" << x0 << ") = " << function(x0) << std::endl;
    out << "Погрешность " << result - function(x0) << std::endl;

    out.close();
}

void task9(const std::string &inputFileName, const std::string &outputFileName)
{

}

void task10(const std::string &inputFileName, const std::string &outputFileName)
{

}

void task11(const std::string &inputFileName, const std::string &outputFileName)
{
    std::ofstream out(outputFileName);

    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 2.0, 3.4142, 4.7321, 6.0};
    double x0 = 2.0;

    int index;
    for (int i = 0; i < x.size(); ++i) if (std::abs(x[i] - x0) < 0.000001) index = i - 1;

    out << "Левосторонняя производная " << "y'(" << x0 << ") = " << Differentiation::firstDerivative(x, y, index) << std::endl;
    out << "Правосторонняя производная " << "y'(" << x0 << ") = " << Differentiation::firstDerivative(x, y, index + 1) << std::endl;
    out << "Первая производная " << "y'(" << x0 << ") = " << Differentiation::firstDerivativeApproximation(x, y, x0, index) << std::endl;
    out << "Вторая производная " << "y''(" << x0 << ") = " << Differentiation::secondDerivative(x, y, index) << std::endl;

    out.close();
}



void task12(const std::string &inputFileName, const std::string &outputFileName)
{
    std::ofstream out(outputFileName);
    out.precision(9);

    double a = 0.0;
    double b = 2.0;

    std::vector<double> h = { 0.5, 0.25 };

    std::vector<std::string> methodNames =
            {
                "Результат метода прямоугольников при h = ",
                "Результат метода трапеций при h = ",
                "Результат метода Симпсона при h = "
            };

    auto function = [](double x) { return (x * x) / (pow(x, 4) + 256); };

    std::vector<AbstractIntegrater*> integrators =
            {
                new RectangleMethod(),
                new TrapezoidMethod(),
                new SimpsonMethod()
            };

    for (int i = 0; i < methodNames.size(); ++i)
    {
        double result_h1 = integrators[i]->integrate(function, a, b, h[0]);
        out << methodNames[i] << h[0] << ": " << result_h1 << std::endl;

        double result_h2 = integrators[i]->integrate(function, a, b, h[1]);
        out << methodNames[i] << h[1] << ": " << result_h2 << std::endl;

        double result = AbstractIntegrater::rungeRoombergMethod(result_h1, result_h2, 2);
        out << "Погрешность: " << std::abs(result - result_h2) << std::endl;

        out << "---------------------------------------------------" << std::endl;
    }

    for (auto x : integrators) delete x;

    out.close();
}

int main() {

    const std::string inputFileName = "../FilesWithResults/input.txt";
    const std::string outputFileName = "../FilesWithResults/outfile.txt";

//    task1(inputFileName, outputFileName);
//    task2(inputFileName, outputFileName);
//    task3(inputFileName, outputFileName);
//    task4(inputFileName, outputFileName);
//    task5(inputFileName, outputFileName);
//    task6(inputFileName, outputFileName);
//    task7(inputFileName, outputFileName);
//    task8(inputFileName, outputFileName);
//    task9(inputFileName, outputFileName);
//    task10(inputFileName, outputFileName);
//    task11(inputFileName, outputFileName);
    task12(inputFileName, outputFileName);

    return 0;
}

