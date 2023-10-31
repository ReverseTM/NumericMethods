#ifndef NUMERICALMETHODS_LU_DECOMPOSITION_H
#define NUMERICALMETHODS_LU_DECOMPOSITION_H

#include "../Matrix/Matrix.h"

class LU_decomposition final
{
    void lu_decomposition(Matrix &A, Matrix& L, Matrix& U, Matrix& P, int &count_per);

    [[nodiscard]] std::vector<double> solve(std::vector<double> b);

    [[nodiscard]] double get_determinant();

    [[nodiscard]] Matrix get_inverse();
};


#endif //NUMERICALMETHODS_LU_DECOMPOSITION_H
