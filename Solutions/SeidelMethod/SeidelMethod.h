#ifndef ZEIDEL_METHOD_H
#define ZEIDEL_METHOD_H

#include "../AbstractSolution.h"

class SeidelMethod final : public AbstractSolution
{

private:

    double epsilon;

    long long iterations = 0;

public:

    SeidelMethod(Matrix &A, double eps);

    std::vector<double> solution(std::vector<double> vector) override;

    long long getCountIterations() const;

private:

    bool isCorrectMatrix() const;


};


#endif //ZEIDEL_METHOD_H
