#include "LagrangePolynomial.h"

std::tuple<std::string, double> LagrangePolynomial::interpolation(
    double (*function)(double),
    const std::vector<double> & x,
    double x0)
{
    int size = x.size();
    double result = 0.0;
    std::string polynomial("L" + std::to_string(size - 1) + "(x) = ");

    for (int i = 0; i < size; ++i)
    {
        double P = 1.0;
        double w = 1.0;

        std::string tmp;

        for (int j = 0; j < size; ++j)
        {
            if (i != j)
            {
                P *= (x0 - x[j]) / (x[i] - x[j]);
                w *= (x[i] - x[j]);

                tmp += "(x" + getNumber(-x[j]) + ")";
            }
        }

        polynomial += getNumber(function(x[i]) / w) + tmp + "\n";
        result += P * function(x[i]);
    }

    return std::make_tuple(polynomial, result);
}

std::string LagrangePolynomial::getNumber(double x)
{
    std::string str = std::to_string(x);
    str = str.substr(0, str.find_last_not_of('0') + 1);

    return str;
}