#ifndef LAGRANGE_POLYNOMIAL_H
#define LAGRANGE_POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <tuple>

class LagrangePolynomial final
{

public:

    static std::tuple<std::string, double> interpolation(
            double (*function)(double),
            const std::vector<double> & x,
            double x0
            );

private:

    static std::string getNumber(double x);

};


#endif //LAGRANGE_POLYNOMIAL_H
