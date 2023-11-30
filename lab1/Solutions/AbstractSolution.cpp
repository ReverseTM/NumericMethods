#include "AbstractSolution.h"

AbstractSolution::AbstractSolution(Matrix &A) :
    matrix(A),
    rows(A.rows),
    cols(A.cols)
{

}
