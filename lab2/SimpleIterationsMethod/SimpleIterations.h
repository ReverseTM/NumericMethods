#ifndef SIMPLE_ITERATIONS_METHOD1_H
#define SIMPLE_ITERATIONS_METHOD1_H

#include <iostream>
#include <tuple>

class SimpleIterationsMethod final
{

public:

    static std::tuple<double, int> solve(
            double (*function)(double),
            double x0,
            double epsilon = 1e-3,
            int maxIterations = 1000);

};


#endif //SIMPLE_ITERATIONS_METHOD1_H
