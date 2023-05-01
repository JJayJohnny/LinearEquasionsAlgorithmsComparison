#include "matrix.h"
#include <stdexcept>

Matrix::Matrix(int rows, int columns, std::vector<std::vector<double>> matrix){
    this->rows = rows;
    this->columns = columns;
    this->matrix = matrix;
}

Matrix::Matrix(int rows, int columns){
    this->rows=rows;
    this->columns=columns;
    this->matrix.resize(rows);
    for(int i=0; i<rows; i++)
        this->matrix[i].resize(columns);
}

void Matrix::Set(int i, int j, double value){
    this->matrix[i][j] = value;
}

double Matrix::Get(int i, int j){
    return this->matrix[i][j];
}

int Matrix::GetRows(){
    return this->rows;
}
int Matrix::GetColumns(){
    return this->columns;
}

Matrix Matrix::Add(Matrix second){
    Matrix result(rows, columns);
    if(second.GetColumns() != columns || second.GetRows() != rows){
        throw std::runtime_error("Matrixes are not compatible for addition");
    }
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            result.Set(i, j, Get(i,j) + second.Get(i, j));
    return result;
}

Matrix Matrix::Subtract(Matrix second){
    Matrix result(rows, columns);
    if(second.GetColumns() != columns || second.GetRows() != rows){
        throw std::runtime_error("Matrixes are not compatible for subtraction");
    }
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            result.Set(i, j, matrix[i][j] - second.Get(i, j));
    return result;
}

Matrix Matrix::Dot(Matrix second){
    Matrix result(rows, second.GetColumns());
    if(columns != second.GetRows())
        throw std::runtime_error("Matrixes are not compatible for multiplication");
    for(int i=0; i<rows; i++){
        for(int j=0; j<second.GetColumns(); j++){
            double sum = 0;
            for(int k=0; k<columns; k++){
                sum+=Get(i, k)*second.Get(k, j);
            }
            result.Set(i, j, sum);
        }
    }
    return result;
}

Matrix Matrix::AddScalar(double scalar){
    Matrix result(rows, columns, matrix);
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            result.Set(i, j, result.Get(i, j)+scalar);
        }
    }
    return result;
}
Matrix Matrix::MultiplyByScalar(double scalar){
    Matrix result(rows, columns, matrix);
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            result.Set(i, j, result.Get(i, j)*scalar);
        }
    }
    return result;
}

std::string Matrix::ToString(){
    std::string s = "";
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            s+=std::to_string(matrix[i][j]);
            s+=" ";
        }
        s+="\n";
    }
    return s;
}

Matrix Matrix::operator+(const Matrix& b){
    return Add(b);
}

Matrix Matrix::operator+(const double& b){
    return AddScalar(b);
}

Matrix Matrix::operator-(const Matrix& b){
    return Subtract(b);
}

Matrix Matrix::operator*(const Matrix& b){
    return Dot(b);
}

Matrix Matrix::operator*(const double& b){
    return MultiplyByScalar(b);
}