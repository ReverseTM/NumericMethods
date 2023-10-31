#ifndef NUMERICALMETHODS_LU_DECOMPOSITION_H
#define NUMERICALMETHODS_LU_DECOMPOSITION_H

#include "../Matrix/Matrix.h"

class LU_decomposition final
{
private:

    Matrix _A;
    Matrix _L;
    Matrix _U;
    Matrix _P;
    int _count_permutations;

public:

    explicit LU_decomposition(Matrix &A);

public:

    static void lu_decomposition(Matrix &A, Matrix& L, Matrix& U, Matrix& P, int &count_per);

    [[nodiscard]] std::vector<double> solve(std::vector<double> b);

    [[nodiscard]] double get_determinant();

    [[nodiscard]] Matrix get_inverse();

    Matrix get_L() const;
    Matrix get_U() const;
    Matrix get_P() const;
    Matrix get_PLU() const;
};


#endif //NUMERICALMETHODS_LU_DECOMPOSITION_H
