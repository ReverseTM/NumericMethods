#ifndef SIMPSON_METHOD_H
#define SIMPSON_METHOD_H

#include "../AbstractIntegrater.h"

class SimpsonMethod final : public AbstractIntegrater
{
    double integrate(double (*function)(double), double a, double b, double h) override;
};


#endif //SIMPSON_METHOD_H
