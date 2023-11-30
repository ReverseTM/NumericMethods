#include "NewtonMethod.h"

std::tuple<double, int> NewtonMethod::solve(
        double (*function)(double),
        double (*derivative)(double),
        double x0,
        double epsilon,
        int maxIterations)
{
    int iteration = 0;

    double previous_x = x0;

    for (; iteration < maxIterations; ++iteration)
    {
        double delta_x = function(previous_x) / derivative(previous_x);
        previous_x = previous_x - delta_x;

        if (std::abs(delta_x) < epsilon) return std::make_tuple(previous_x, iteration);
    }

    throw std::runtime_error("Mat iterations");
}