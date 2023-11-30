#include "SISystem.h"

std::tuple<std::vector<double>, int> SISystem::solve(
        std::vector<double (*)(std::vector<double>)> functions,
        double x0,
        double epsilon,
        int maxIterations)
{
    std::vector<double> x(functions.size(), x0);

    int iteration = 0;

    while (iteration < maxIterations)
    {
        std::vector<double> xNext(functions.size());

        for (int i = 0; i < functions.size(); ++i)
        {
            xNext[i] = functions[i](x);
        }


        double max = 0;
        for (int i = 0; i < xNext.size(); ++i)
        {
            max = std::max(max, std::abs(xNext[i] - x[i]));
        }

        if (max < epsilon)
        {
            return std::make_tuple(xNext, iteration);
        }

        x = xNext;
        iteration++;
    }

    throw std::runtime_error("Max iterations");
}