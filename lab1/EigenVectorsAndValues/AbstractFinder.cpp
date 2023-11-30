#include "AbstractFinder.h"

AbstractFinder::AbstractFinder(Matrix &A) :
    matrix(A),
    rows(A.rows),
    cols(A.cols)
{

}