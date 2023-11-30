#ifndef NMSYSTEM_H
#define NMSYSTEM_H

#include <iostream>
#include <vector>
#include <tuple>
#include "../../lab1/Matrix/Matrix.h"
#include "../../lab1/Solutions/LU_decomposition/LUDecomposition.h"

class NMSystem final
{

public:

    static std::tuple<std::vector<double>, int> solve(
            std::vector<double (*)(std::vector<double>)> functions,
            std::vector<double (*)(std::vector<double>)> derivatives,
            double x0,
            double epsilon = 1e-3,
            int maxIterations = 1000);

};


#endif //NMSYSTEM_H
