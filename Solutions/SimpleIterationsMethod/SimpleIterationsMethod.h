#ifndef SIMPLE_ITERATIONS_METHOD_H
#define SIMPLE_ITERATIONS_METHOD_H

#include "../AbstractSolution.h"

class SimpleIterationsMethod final : public AbstractSolution
{

private:

    long long iterations = 0;

public:

    explicit SimpleIterationsMethod(Matrix &A);

public:

    std::vector<double> solution(std::vector<double> b, double epsilon = 0.1, int maxIterations = 1000) override;

public:

    [[nodiscard]] long long getCountIterations() const;

private:

    bool isCorrectMatrix() const;

};


#endif //SIMPLE_ITERATIONS_METHOD_H
