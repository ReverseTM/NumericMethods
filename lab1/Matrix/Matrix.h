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
    friend class SimpleIterations;
    friend class SeidelMethod;
    friend class AbstractFinder;
    friend class RotationMethod;
    friend class QRDecomposition;
    friend class NMSystem;

private:

    std::vector<std::vector<double>> data;

    int rows;
    int cols;

public:

    explicit Matrix();

    explicit Matrix(int, int);

    void inputMatrix();

    void inputMatrixFromFile(std::ifstream & filename);

    int getRows() const;

    int getCols() const;

    std::vector<double> getRow(int rowIndex) const;

    std::vector<double> getCol(int colIndex) const;

    void printRow(std::ostream &out, int rowIndex) const;

    void printCol(std::ostream &out,int colIndex) const;

public:

    friend std::ostream& operator<<(std::ostream &out, const Matrix &matrix);

    friend Matrix operator*(const Matrix &A, const Matrix &B);

    friend Matrix operator*(const Matrix &A, double x);

    friend Matrix operator*(double x, const Matrix &A);

    void operator*() const;

    friend std::vector<double> operator*(const Matrix &A, const std::vector<double> &vector);

    friend std::vector<double> operator*(const std::vector<double> &vector, const Matrix &B);

    friend Matrix operator-(const Matrix &A, const Matrix &B);

    friend Matrix operator+(const Matrix &A, const Matrix &B);

    friend bool operator==(const Matrix &A, const Matrix &B);

    friend bool operator!=(const Matrix &A, const Matrix &B);

    friend Matrix operator~(const Matrix &A);

private:

    static void transform_l(Matrix &L, int i, int j);

    static int find_pivot_row(const Matrix & matrix, int col);

    static void swap_row(Matrix & matrix, int i, int j);
};

#endif //MATRIX_H