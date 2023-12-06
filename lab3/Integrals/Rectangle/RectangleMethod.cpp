#include "RectangleMethod.h"

double RectangleMethod::integrate(
    double (*function)(double),
    double a,
    double b,
    double h)
{
    int n = (b - a) / h;

    double sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double x = a + (i + 0.5) * h;
        sum += function(x);
    }

    return sum * h;
}