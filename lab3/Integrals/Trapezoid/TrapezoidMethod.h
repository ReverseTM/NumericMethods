#ifndef TRAPEZOID_METHOD_H
#define TRAPEZOID_METHOD_H

#include "../AbstractIntegrater.h"

class TrapezoidMethod final : public AbstractIntegrater
{
    double integrate(double (*function)(double), double a, double b, double h) override;
};


#endif //TRAPEZOID_METHOD_H
