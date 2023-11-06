#include "RotationMethod.h"

RotationMethod::RotationMethod(Matrix &A) :
    AbstractFinder(A)
{

}

std::tuple<std::vector<double>, Matrix> RotationMethod::find(Matrix &A, double epsilon, int maxIterations)
{

}