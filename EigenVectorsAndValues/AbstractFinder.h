#ifndef ABSTRACT_FINDER_H
#define ABSTRACT_FINDER_H

#include <iostream>
#include <vector>
#include <tuple>

#include "../Matrix/Matrix.h"

class AbstractFinder
{

protected:

    Matrix matrix;

    int rows;
    int cols;

public:

    explicit AbstractFinder(Matrix &A);

    virtual ~AbstractFinder() = default;

public:

    virtual std::tuple<std::vector<double>, Matrix> find(Matrix &A, double epsilon, int maxIterations) = 0;

};


#endif //ABSTRACT_FINDER_H
