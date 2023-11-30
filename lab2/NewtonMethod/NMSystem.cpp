#include "NMSystem.h"

std::tuple<std::vector<double>, int> NMSystem::solve(
    std::vector<double (*)(std::vector<double>)> functions,
    std::vector<double (*)(std::vector<double>)> derivatives,
    double x0,
    double epsilon,
    int maxIterations)
{
    std::vector<double> x(functions.size(), x0);
    std::vector<double> xNext(functions.size());

    int iteration = 0;

    while (iteration < maxIterations)
    {
        Matrix jacobian(functions.size(), functions.size());
        for (int i = 0; i < jacobian.rows; ++i)
        {
            for (int j = 0; j < jacobian.cols; ++j)
            {
                jacobian.data[i][j] = derivatives[functions.size() * i + j](x);
            }
        }

        std::vector<double> b(functions.size());
        for (int i = 0; i < b.size(); ++i) b[i] = -functions[i](x);

        try
        {
            LUDecomposition solver(jacobian);
            auto delta_x = solver.solution(b);

            for (int i = 0; i < xNext.size(); ++i) xNext[i] = x[i] + delta_x[i];

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

        }
        catch (std::exception &ex)
        {
            throw ex;
        }

        iteration++;
    }

    throw std::runtime_error("Max iterations");
}