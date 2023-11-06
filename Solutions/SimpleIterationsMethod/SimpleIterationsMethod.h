#ifndef SIMPLE_ITERATIONS_METHOD_H
#define SIMPLE_ITERATIONS_METHOD_H

#include "../AbstractSolution.h"

class SimpleIterationsMethod final : public AbstractSolution
{

private:

    double epsilon;

    long long iterations = 0;

public:

    explicit SimpleIterationsMethod(Matrix &A, double eps);

public:

    std::vector<double> solution(std::vector<double> vector) override;

public:

    [[nodiscard]] long long getCountIterations() const;

private:

    bool isCorrectMatrix() const;

};


#endif //SIMPLE_ITERATIONS_METHOD_H
