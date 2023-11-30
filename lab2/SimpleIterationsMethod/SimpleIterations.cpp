#include "SimpleIterations.h"

std::tuple<double, int> SimpleIterationsMethod::solve(
        double (*function)(double),
        double x0,
        double epsilon,
        int maxIterations)
{
    int iteration = 0;

    double previous_x = x0;

    for (; iteration < maxIterations; ++iteration)
    {
        double next_x = function(previous_x);

        if (std::abs(next_x - previous_x) < epsilon) return std::make_tuple(next_x, iteration);

        previous_x = next_x;
    }

    throw std::runtime_error("Max iterations!");
}