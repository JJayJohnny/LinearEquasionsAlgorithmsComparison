#include <iostream>
#include <vector>
#include <chrono>
#include <matplotlibcpp.h>
#include "matrix.h"

#define PLOT_DIR "plots/"
#define N 567 //9*7*9
#define A1 11
#define A2 -1
#define A3 -1
#define F 8
//index 188679

namespace plt = matplotlibcpp;

Matrix CalculateJacobi(Matrix& A, Matrix& b, double accuracy, std::vector<double>& iterations, std::vector<double>& residuum){
    Matrix L = A.GetTriangleL();
    Matrix D = A.GetDiagonal();
    Matrix U = A.GetTriangleU();
    Matrix r(A.GetRows(), 1);
    r.FillOnes();
    double error = (A*r-b).CalculateNorm();
    double i=0;
    //these dont change through the algorithm so I preallocated them for speed
    Matrix rParam = (D*(-1)).InverseDiagonalMatrix() * (L+U);
    Matrix freeParam = D.InverseDiagonalMatrix()*b;
    while(error > accuracy){
        r = rParam*r + freeParam;
        error = (A*r-b).CalculateNorm();
        i++;
        iterations.push_back(i);
        residuum.push_back(error);
    }
    return r;
}

Matrix CalculateGauss(Matrix& A, Matrix& b, double accuracy, std::vector<double>& iterations, std::vector<double>& residuum){
    Matrix L = A.GetTriangleL();
    Matrix D = A.GetDiagonal();
    Matrix U = A.GetTriangleU();
    Matrix r(A.GetRows(), 1);
    r.FillOnes();
    double error = (A*r-b).CalculateNorm();
    double i=0;
    //these dont change through the algorithm so I preallocated them for speed
    Matrix rParam = (D+L)*(-1);
    Matrix freeParam = (D+L).ForwardSubstitution(b);
    while(error > accuracy){
        r=rParam.ForwardSubstitution(U*r) + freeParam;
        error = (A*r-b).CalculateNorm();
        i++;
        iterations.push_back(i);
        residuum.push_back(error);
    }
    return r;
}

void Example1(){
    std::cout<<"Zadanie B:"<<"\n";
    Matrix A(N, N);
    A.CreateDiagonal(A1, A2, A3);
    Matrix b(N, 1);
    for(int i=0; i<N; i++){
        b.Set(i, 0, sin(i*(F+1)));
    }
    std::vector<double> residuumJacobi;
    std::vector<double> iterationsJacobi;
    auto startJacobi = std::chrono::high_resolution_clock::now();
    Matrix rJacobi = CalculateJacobi(A, b, 1e-9, iterationsJacobi, residuumJacobi);
    auto stopJacobi = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Jacobiego-"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterationsJacobi.size()<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopJacobi-startJacobi).count()<<"ms"<<"\n";
    //metoda Gaussa-Seidla
    std::vector<double> residuumGauss;
    std::vector<double> iterationsGauss;
    auto startGauss = std::chrono::high_resolution_clock::now();
    Matrix rGauss = CalculateGauss(A, b, 1e-9, iterationsGauss, residuumGauss);
    auto stopGauss = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Gaussa-Seidla-"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterationsGauss.size()<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopGauss-startGauss).count()<<"ms"<<"\n";
    //generowanie wykresu
    plt::figure();
    plt::figure_size(700, 500);
    plt::named_semilogy("Metoda Jacobiego", iterationsJacobi, residuumJacobi);
    plt::named_semilogy("Metoda Gaussa-Seidla", iterationsGauss, residuumGauss);
    plt::title("Porownanie bledu rezydualnego w kolejnych iteracjach");
    plt::xlabel("Iteracja");
    plt::ylabel("Norma bledu rezydualnego");
    plt::legend();
    plt::save("plots/ZadanieBResiduum.png");
}

int main(){
    Matrix m(5, 2, {{5.0, 1.0},{-3.0, 3.0},{4.0, -1.0},{-1.0, 2.0},{0.0, 7.0}});
    Matrix d(2, 7, {{1.0, 7.0, 0.0, -1.0, 9.0, 0.0, 5.0},{2.0, -6.0, -1.0, 10.0, 2.0, -3.0, 1.0}});
    

    //zadanie A i B
    Example1();
    
    return 0;
}