#ifndef NUMERICALMETHODS_SWEEP_METHOD_H
#define NUMERICALMETHODS_SWEEP_METHOD_H

#include "../Matrix/Matrix.h"

class sweep_method final
{
private:

    Matrix matrix;

    std::vector<double> P;
    std::vector<double> Q;

    int rows;
    int cols;

public:

    sweep_method(Matrix &A);

    std::vector<double> solution(std::vector<double> b);

private:

    bool is_correct_matrix() const;

};


#endif //NUMERICALMETHODS_SWEEP_METHOD_H
