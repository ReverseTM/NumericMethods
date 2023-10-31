#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class Matrix final {

    friend class LU_decomposition;

private:

    std::vector<std::vector<double>> data;

    int rows;
    int cols;

    std::string debug_file;
    std::string output_file;

public:

    Matrix(int, int);

    void inputMatrix();

    void inputMatrixFromFile(std::ifstream & filename);

    void lu_decomposition(Matrix& L, Matrix& U, Matrix& P, int &count_per) const;

    [[nodiscard]] std::vector<double> solve(std::vector<double> b) const;

    [[nodiscard]] double get_determinant() const;

    [[nodiscard]] Matrix get_inverse() const;

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
