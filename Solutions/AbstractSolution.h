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

    explicit AbstractSolution(Matrix &A);

    virtual ~AbstractSolution() = default;

public:

    virtual std::vector<double> solution(std::vector<double> vector, double epsilon = 0.1, int maxIterations = 1000) = 0;
};

#endif //ABSTRACT_SOLUTION_H
