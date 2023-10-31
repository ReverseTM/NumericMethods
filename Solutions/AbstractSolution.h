#ifndef ABSTRACT_SOLUTION_H
#define ABSTRACT_SOLUTION_H
#include <iostream>
#include <vector>

#include "../Matrix/Matrix.h"

class AbstractSolution
{

protected:

    Matrix matrix;

    int rows;
    int cols;

public:

    virtual std::vector<double> solution(std::vector<double> vector) = 0;
    virtual ~AbstractSolution() = default;
};

#endif //ABSTRACT_SOLUTION_H
