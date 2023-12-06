#include "Differentiation.h"

double Differentiation::firstDerivative(
        const std::vector<double> &x,
        const std::vector<double> &y,
        int i)
{
    return (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

double Differentiation::firstDerivativeApproximation(
    const std::vector<double> &x,
    const std::vector<double> &y,
    double x0,
    int i)
{
    double leftDerivative = firstDerivative(x, y, i);
    double rightDerivative = firstDerivative(x, y, i + 1);

    return leftDerivative + (rightDerivative - leftDerivative) / (x[i + 2] - x[i]) * (2 * x0 - x[i] - x[i + 1]);
}

double Differentiation::secondDerivative(const std::vector<double> &x, const std::vector<double> &y, int i)
{
    double leftDerivative = firstDerivative(x, y, i);
    double rightDerivative = firstDerivative(x, y, i + 1);

    return 2 * ((rightDerivative - leftDerivative) / (x[i + 2] - x[i]));
}