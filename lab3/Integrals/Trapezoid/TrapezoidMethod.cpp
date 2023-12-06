#include "TrapezoidMethod.h"

double TrapezoidMethod::integrate(
    double (*function)(double),
    double a,
    double b,
    double h)
{
    int n = (b - a) / h;

    double sum = (function(a) + function(b)) / 2;

    for (int i = 1; i < n; ++i)
    {
        double x = a + i * h;
        sum += function(x);
    }

    return sum * h;
}