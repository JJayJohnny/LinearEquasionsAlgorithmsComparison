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

Matrix CalculateJacobi(Matrix& A, Matrix& b, double accuracy, std::vector<int>& iterations, std::vector<double>& residuum){
    Matrix L = A.GetTriangleL();
    Matrix D = A.GetDiagonal();
    Matrix U = A.GetTriangleU();
    Matrix r(A.GetRows(), 1);
    r.FillOnes();
    double error = (A*r-b).CalculateNorm();
    int i=0;
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

int main(){
    Matrix m(5, 2, {{5.0, 1.0},{-3.0, 3.0},{4.0, -1.0},{-1.0, 2.0},{0.0, 7.0}});
    Matrix d(2, 7, {{1.0, 7.0, 0.0, -1.0, 9.0, 0.0, 5.0},{2.0, -6.0, -1.0, 10.0, 2.0, -3.0, 1.0}});
    

    //zadanie A
    Matrix A(N, N);
    A.CreateDiagonal(A1, A2, A3);
    Matrix b(N, 1);
    for(int i=0; i<N; i++){
        b.Set(i, 0, sin(i*(F+1)));
    }
    std::vector<double> residuum;
    std::vector<int> iterations;
    auto start = std::chrono::high_resolution_clock::now();
    Matrix r = CalculateJacobi(A, b, 1e-9, iterations, residuum);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout<<"Metoda Jacobiego"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterations.size()<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count()<<"\n";
    plt::semilogy(iterations, residuum);
    plt::title("Wykres bledu rezydualnego dla algorytmu Jacobiego");
    plt::xlabel("Iteracja");
    plt::ylabel("Norma bledu rezydualnego");
    plt::save("plots/JacobiResiduum.png");
    return 0;
}