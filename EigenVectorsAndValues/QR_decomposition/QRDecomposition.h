#ifndef QR_DECOMPOSITION_H
#define QR_DECOMPOSITION_H

#include "../AbstractFinder.h"

class QRDecomposition final
{

private:

    Matrix matrix;
    Matrix Q;
    Matrix R;

    int rows;
    int cols;

    long long iterations = 0;

public:

    explicit QRDecomposition(Matrix &A);

public:

    std::vector<std::string> find(double epsilon, int maxIterations = 1000);

    [[nodiscard]] long long getCountIterations() const;

    Matrix get_R() const;
    Matrix get_Q() const;

private:

    void qr_decomposition();

    static bool eigenValuesReady(const Matrix &A, double epsilon);

};

#endif //QR_DECOMPOSITION_H
