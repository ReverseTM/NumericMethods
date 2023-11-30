#ifndef SISYSTEM_H
#define SISYSTEM_H

#include <iostream>
#include <vector>
#include <tuple>

class SISystem final
{

public:

    static std::tuple<std::vector<double>, int> solve(
            std::vector<double (*)(std::vector<double>)> functions,
            double x0,
            double epsilon = 1e-3,
            int maxIterations = 1000);

};


#endif //SISYSTEM_H
