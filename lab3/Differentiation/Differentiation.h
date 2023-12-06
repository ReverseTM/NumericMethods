#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H

#include <iostream>
#include <vector>
#include <tuple>

class Differentiation final
{

public:

    static double firstDerivative(const std::vector<double> & x, const std::vector<double> & y, int i);

    static double firstDerivativeApproximation(const std::vector<double> & x, const std::vector<double> & y, double x0, int i);

    static double secondDerivative(const std::vector<double> & x, const std::vector<double> & y, int i);

};


#endif //DIFFERENTIATION_H
