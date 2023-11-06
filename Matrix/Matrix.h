#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class Matrix final {

    friend class AbstractSolution;
    friend class LUDecomposition;
    friend class SweepMethod;
    friend class SimpleIterationsMethod;
    friend class SeidelMethod;
    friend class AbstractFinder;
    friend class RotationMethod;

private:

    std::vector<std::vector<double>> data;

    int rows;
    int cols;

public:

    explicit Matrix();

    explicit Matrix(int, int);

    void inputMatrix();

    void inputMatrixFromFile(std::ifstream & filename);

public:

    friend std::ostream& operator<<(std::ostream &out, const Matrix &matrix);

    friend Matrix operator*(const Matrix &A, const Matrix &B);

    friend Matrix operator*(const Matrix &A, int x);

    friend Matrix operator*(int x, const Matrix &A);

    void operator*() const;

    friend std::vector<double> operator*(const Matrix &A, const std::vector<double> &vector);

    friend std::vector<double> operator*(const std::vector<double> &vector, const Matrix &B);

    friend Matrix operator-(const Matrix &A, const Matrix &B);

    friend Matrix operator+(const Matrix &A, const Matrix &B);

    friend Matrix operator~(const Matrix &A);

private:

    static void transform_l(Matrix &L, int i, int j);

    static int find_pivot_row(const Matrix & matrix, int col);

    static void swap_row(Matrix & matrix, int i, int j);
};

#endif //MATRIX_H