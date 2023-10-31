#ifndef SWEEP_METHOD_H
#define SWEEP_METHOD_H

#include "../AbstractSolution.h"

class SweepMethod final : public AbstractSolution
{

private:

    std::vector<double> P;
    std::vector<double> Q;

public:

    explicit SweepMethod(Matrix &A);

    std::vector<double> solution(std::vector<double> b) override;

private:

    [[nodiscard]] bool is_correct_matrix() const;

    [[nodiscard]] std::vector<double> get_P() const;
    [[nodiscard]] std::vector<double> get_Q() const;


};


#endif //SWEEP_METHOD_H
