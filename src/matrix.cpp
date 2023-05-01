#include "matrix.h"

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

Matrix* Matrix::Add(Matrix second){
    Matrix* result = nullptr;
    if(second.GetColumns() != columns || second.GetRows() != rows){
        return result;
    }
    result = new Matrix(rows, columns);
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            result->Set(i, j, Get(i,j) + second.Get(i, j));
    return result;
}

Matrix* Matrix::Substract(Matrix second){
    Matrix* result = nullptr;
    if(second.GetColumns() != columns || second.GetRows() != rows){
        return result;
    }
    result = new Matrix(rows, columns);
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            result->Set(i, j, matrix[i][j] - second.Get(i, j));
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