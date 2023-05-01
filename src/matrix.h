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
        void Set(int i, int j, double value);
        double Get(int i, int j);
        int GetRows();
        int GetColumns();
        Matrix Add(Matrix second);
        Matrix Subtract(Matrix second);
        Matrix Dot(Matrix second);
        std::string ToString();
};