#include "LU_decomposition.h"

void LU_decomposition::lu_decomposition(Matrix &A, Matrix &L, Matrix &U, Matrix &P, int &count)
{
    if (L.rows != L.cols)
    {
        return;
    }

    U = *this;

    for (int i = 0; i < U.rows; i++)
    {
        L.data[i][i] = 1.0;
        P.data[i][i] = 1.0;
    }

    for (int i = 0; i < U.rows - 1; i++)
    {
        int pivot = Matrix::find_pivot_row(U, i);
        if (pivot != i)
        {
            Matrix::swap_row(U, i, pivot);
            Matrix::swap_row(P, i, pivot);
            Matrix::transform_l(L, i, pivot);
            count++;
        }

        for (int j = i + 1; j < U.rows; j++)
        {
            double factor = U.data[j][i] / U.data[i][i];
            L.data[j][i] = factor;

            for (int k = i; k < U.rows; k++)
            {
                U.data[j][k] -= U.data[i][k] * factor;
            }
        }
    }

    P = ~P;
}
