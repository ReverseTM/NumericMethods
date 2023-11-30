#ifndef NEWTON_METHOD_H
#define NEWTON_METHOD_H

#include <iostream>
#include <tuple>

class NewtonMethod final
{

public:

    static std::tuple<double, int> solve(
            double (*function)(double),
            double (*derivative)(double),
            double x0,
            double epsilon,
            int maxIterations = 1000);

};


#endif //NEWTON_METHOD_H
