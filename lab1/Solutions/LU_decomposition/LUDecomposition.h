#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include "../AbstractSolution.h"

class LUDecomposition final : public AbstractSolution
{

private:

    Matrix L;
    Matrix U;
    Matrix P;

    int count_permutations = 0;

public:

    explicit LUDecomposition(Matrix &A);

public:

    [[nodiscard]] std::vector<double> solution(std::vector<double> vector, double epsilon = 0.1, int maxIterations = 1000) override;

    [[nodiscard]] double get_determinant();

    [[nodiscard]] Matrix get_inverse();

    [[nodiscard]] Matrix get_L() const;
    [[nodiscard]] Matrix get_U() const;
    [[nodiscard]] Matrix get_P() const;
    [[nodiscard]] Matrix get_PLU() const;

private:

    void lu_decomposition();

};

#endif //LU_DECOMPOSITION_H
