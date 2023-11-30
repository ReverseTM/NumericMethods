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

    std::tuple<Matrix, std::vector<double>> find(double epsilon, int maxIterations) override;

    [[nodiscard]] long long getCountIterations() const;

private:

    bool isCorrectMatrix() const;

};

#endif //ROTATION_METHOD_H
