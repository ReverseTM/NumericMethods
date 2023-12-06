#include "NewtonPolynomial.h"

std::tuple<std::string, double> NewtonPolynomial::interpolation(
    double (*function)(double),
    const std::vector<double> & x,
    double x0)
{
    int size = x.size();

    std::vector<double> y(size);
    for (int i = 0; i < size; ++i) y[i] = function(x[i]);

    double result = y[0];
    std::string polynomial("P" + std::to_string(size - 1) + "(x) = " + getNumber(result) + "\n");

    for (int i = 1; i < size; ++i)
    {
        double difference = dividedDifference(x, y, i);

        double mult = 1.0;
        std::string tmp;

        for (int j = 0; j < i; ++j)
        {
            mult *= (x0 - x[j]);

            tmp += "(x" + getNumber(-x[j]) + ")";

        }

        if (difference >= 0) polynomial += "+";

        polynomial += getNumber(difference) + tmp + "\n";
        result += difference * mult;
    }

    return std::make_tuple(polynomial, result);
}

double NewtonPolynomial::dividedDifference(const std::vector<double> & x, const std::vector<double> &y, int k)
{
    double result = 0.0;

    for (int i = 0; i <= k; ++i)
    {
        double mult = 1.0;
        for (int j = 0; j <= k; ++j)
        {
            if (j != i) mult *= x[i] - x[j];
        }

        result += y[i] / mult;
    }

    return result;
}

std::string NewtonPolynomial::getNumber(double x)
{
    std::string str = std::to_string(x);
    str = str.substr(0, str.find_last_not_of('0') + 1);

    return str;
}