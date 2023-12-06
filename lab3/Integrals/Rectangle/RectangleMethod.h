#ifndef RECTANGLE_METHOD_H
#define RECTANGLE_METHOD_H

#include "../AbstractIntegrater.h"

class RectangleMethod final : public AbstractIntegrater
{

public:

    double integrate(double (*function)(double), double a, double b, double h) override;

};


#endif //RECTANGLE_METHOD_H
