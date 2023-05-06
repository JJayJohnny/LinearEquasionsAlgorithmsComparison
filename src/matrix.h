#pragma once
#include <vector>
#include <string>

class Matrix{
    private:
        int rows, columns;
        std::vector<std::vector<double>> matrix;
    public:
        Matrix(int rows, int columns, std::vector<std::vector<double>>);
        Matrix(int rows, int columns);
        Matrix(Matrix& b);
        void Set(int i, int j, double value);
        double Get(int i, int j) const;
        int GetRows() const;
        int GetColumns() const;
        Matrix Add(const Matrix& b);
        Matrix Subtract(const Matrix& b);
        Matrix Dot(const Matrix& b);
        Matrix AddScalar(const double scalar);
        Matrix MultiplyByScalar(const double scalar);
        Matrix Transpose();
        std::string ToString();
        Matrix operator+(const Matrix& b);
        Matrix operator+(const double& b);
        Matrix operator-(const Matrix& b);
        Matrix operator*(const Matrix& b);
        Matrix operator*(const double& b);

        void FillZeros();
        void FillOnes();
        Matrix GetTriangleL() const;
        Matrix GetTriangleU() const;
        Matrix GetDiagonal() const;
        double CalculateNorm();
        void CreateDiagonal(double a1, double a2, double a3);
        Matrix ForwardSubstitution(const Matrix& b);
        Matrix BackwardSubstitution(const Matrix& b);
};