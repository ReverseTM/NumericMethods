#ifndef ZEIDEL_METHOD_H
#define ZEIDEL_METHOD_H

#include "../AbstractSolution.h"

class SeidelMethod final : public AbstractSolution
{

private:

    long long iterations = 0;

public:

    explicit SeidelMethod(Matrix &A);

    std::vector<double> solution(std::vector<double> b, double epsilon = 0.1, int maxIterations = 1000) override;

    long long getCountIterations() const;

private:

    bool isCorrectMatrix() const;


};


#endif //ZEIDEL_METHOD_H
