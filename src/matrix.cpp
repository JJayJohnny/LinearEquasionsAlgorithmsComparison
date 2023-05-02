#include "matrix.h"
#include <stdexcept>
#include <cmath>

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

Matrix::Matrix(Matrix& b){
    this->rows = b.rows;
    this->columns = b.columns;
    this->matrix = b.matrix;
}

void Matrix::Set(int i, int j, double value){
    this->matrix[i][j] = value;
}

double Matrix::Get(int i, int j) const{
    return this->matrix[i][j];
}

int Matrix::GetRows() const{
    return this->rows;
}
int Matrix::GetColumns() const{
    return this->columns;
}

Matrix Matrix::Add(const Matrix& b){
    Matrix result(rows, columns);
    if(b.GetColumns() != columns || b.GetRows() != rows){
        throw std::runtime_error("Matrixes are not compatible for addition");
    }
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            result.Set(i, j, Get(i,j) + b.Get(i, j));
    return result;
}

Matrix Matrix::Subtract(const Matrix& b){
    Matrix result(rows, columns);
    if(b.GetColumns() != columns || b.GetRows() != rows){
        throw std::runtime_error("Matrixes are not compatible for subtraction");
    }
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            result.Set(i, j, matrix[i][j] - b.Get(i, j));
    return result;
}

Matrix Matrix::Dot(const Matrix& b){
    Matrix result(rows, b.GetColumns());
    if(columns != b.GetRows())
        throw std::runtime_error("Matrixes are not compatible for multiplication");
    for(int i=0; i<rows; i++){
        for(int j=0; j<b.GetColumns(); j++){
            double sum = 0;
            for(int k=0; k<columns; k++){
                sum+=Get(i, k)*b.Get(k, j);
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

void Matrix::FillZeros(){
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            Set(i, j, 0);
}

void Matrix::FillOnes(){
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
            Set(i, j, 1);
}

void Matrix::CreateDiagonal(double a1, double a2, double a3){
    FillZeros();
    for(int i=0; i<rows; i++){
        Set(i, i, a1);
        if(i>=1)
            Set(i, i-1, a2);
        if(i>=2)
            Set(i, i-2, a3);
        if(i<columns-1)
            Set(i, i+1, a2);
        if(i<columns-2)
            Set(i, i+2, a3);
    }
}

Matrix Matrix::ForwardSubstitution(const Matrix& b){
    Matrix x(rows, 1);
    x.FillZeros();
    if(rows != b.GetRows())
        throw std::runtime_error("Matrixes are not compatible for forward substitution");
    for(int i=0; i<rows; i++){
        double temp = b.Get(i, 0);
        for(int j=0; j<i; j++){
            temp -= Get(i, j) * x.Get(j, 0);
        }
        x.Set(i, 0, temp / Get(i, i));
    }
    return x;
}

Matrix Matrix::BackwardSubstitution(const Matrix& b){
    Matrix x(rows, 1);
    x.FillZeros();
    if(rows != b.GetRows())
        throw std::runtime_error("Matrixes are not compatible for backward substitution");
    for(int i=rows-1; i>=0; i--){
        double temp = b.Get(i, 0);
        for(int j=i+1; j<columns; j++){
            temp -= Get(i, j) * x.Get(j, 0);
        }
        x.Set(i, 0, temp/Get(i,i));
    }
    return x;
}

Matrix Matrix::GetTriangleL() const{
    Matrix L(rows, columns);
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            if(j<i)
                L.Set(i, j, Get(i, j));
            else
                L.Set(i, j, 0);
        }
    }
    return L;
}
Matrix Matrix::GetTriangleU() const{
    Matrix U(rows, columns);
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            if(j>i)
                U.Set(i, j, Get(i, j));
            else
                U.Set(i, j, 0);
        }
    }
    return U;
}
Matrix Matrix::GetDiagonal() const{
    Matrix D(rows, columns);
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            if(j==i)
                D.Set(i, j, Get(i, j));
            else
                D.Set(i, j, 0);
        }
    }
    return D;

}

Matrix Matrix::Transpose(){
    Matrix result(columns, rows);
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            result.Set(j, i, Get(i, j));
        }
    }
    return result;
}

double Matrix::CalculateNorm(){
    if(rows != 1 && columns != 1)
        throw std::runtime_error("Only vectors can have norm calculated");
    return sqrt(Transpose().Dot(*this).Get(0, 0));
}

Matrix Matrix::InverseDiagonalMatrix(){
    Matrix result(rows, columns, matrix);
    for(int i=0; i<rows; i++){
        result.Set(i, i, 1/Get(i, i));
    }
    return result;
}