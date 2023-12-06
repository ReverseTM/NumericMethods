#include <cmath>
#include "AbstractIntegrater.h"

double AbstractIntegrater::rungeRoombergMethod(
    double result_h1,
    double result_h2,
    int p)
{
    return result_h2 + (result_h2 - result_h1) / (pow(2, p) - 1);
}