#ifndef NEWTON_POLYNOMIAL_H
#define NEWTON_POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <tuple>

class NewtonPolynomial final
{

public:

    static std::tuple<std::string, double> interpolation(
            double (*function)(double),
            const std::vector<double> & x,
            double x0
    );

private:

    static double dividedDifference(const std::vector<double> & x, const std::vector<double> &y, int k);

    static std::string getNumber(double x);
};


#endif //NEWTON_POLYNOMIAL_H
