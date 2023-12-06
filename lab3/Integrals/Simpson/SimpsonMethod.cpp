#include "SimpsonMethod.h"

double SimpsonMethod::integrate(
    double (*function)(double),
    double a,
    double b,
    double h)
{
    int n = (b - a) / h;

    double sum = function(a) + function(b);

    for (int i = 1; i < n; ++i)
    {
        double x = function(a + i * h);
        if (i & 1) sum += 4 * x;
        else sum += 2 * x;
    }

    return sum * h / 3;
}