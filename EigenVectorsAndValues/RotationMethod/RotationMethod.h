#ifndef ROTATION_METHOD_H
#define ROTATION_METHOD_H

#include "../AbstractFinder.h"

class RotationMethod final : public AbstractFinder
{

private:

    long long iterations = 0;

public:

    explicit RotationMethod(Matrix &A);

public:

    std::tuple<std::vector<double>, Matrix> find(Matrix &A, double epsilon, int maxIterations) override;
};

#endif //ROTATION_METHOD_H
