#include <iostream>
#include <vector>
#include <chrono>
#include <matplotlibcpp.h>
#include "matrix.h"

#define PLOT_DIR "plots/"
#define N 979
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
    iterations.push_back(i);
    residuum.push_back(error);
    //these dont change through the algorithm so I preallocated them for speed
    Matrix freeParam = D.ForwardSubstitution(b);
    while(error > accuracy){
        r = D.ForwardSubstitution(((L+U)*r)*(-1)) + freeParam;
        error = (A*r-b).CalculateNorm();
        i++;
        iterations.push_back(i);
        residuum.push_back(error);
        // if(i >= 1000)
        //     break;
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
    iterations.push_back(i);
    residuum.push_back(error);
    //these dont change through the algorithm so I preallocated them for speed
    Matrix rParam = (D+L)*(-1);
    Matrix freeParam = (D+L).ForwardSubstitution(b);
    while(error > accuracy){
        r=rParam.ForwardSubstitution(U*r) + freeParam;
        error = (A*r-b).CalculateNorm();
        i++;
        iterations.push_back(i);
        residuum.push_back(error);
        // if(i >= 1000)
        //     break;
    }
    return r;
}

Matrix LUFactorization(Matrix& A, Matrix& b){
    Matrix U(A);
    Matrix L(A.GetRows(), A.GetColumns());
    L.CreateDiagonal(1, 0, 0);
    for(int k=0; k<A.GetRows()-1; k++){
        for(int j=k+1; j<A.GetColumns(); j++){
            L.Set(j, k, U.Get(j, k)/U.Get(k, k));
            for(int i=k; i<U.GetColumns(); i++){
                U.Set(j, i, U.Get(j, i)-L.Get(j,k)*U.Get(k, i));
            }
        }
    }
    Matrix y = L.ForwardSubstitution(b);
    Matrix x = U.BackwardSubstitution(y);
    return x;
}

void Example1(){
    std::cout<<"\nZadanie B:"<<"\n";
    std::cout<<"N="<<N<<"\n";
    Matrix A(N, N);
    A.CreateDiagonal(A1, A2, A3);
    Matrix b(N, 1);
    for(int i=0; i<N; i++){
        b.Set(i, 0, sin(i*(F+1)));
    }
    //metoda Jacobiego
    std::vector<double> residuumJacobi;
    std::vector<double> iterationsJacobi;
    auto startJacobi = std::chrono::high_resolution_clock::now();
    Matrix rJacobi = CalculateJacobi(A, b, 1e-9, iterationsJacobi, residuumJacobi);
    auto stopJacobi = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Jacobiego-"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterationsJacobi.size()-1<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopJacobi-startJacobi).count()<<"ms"<<"\n";
    std::cout<<"Norma bledu rezydualnego: "<<residuumJacobi[residuumJacobi.size()-1]<<"\n";
    //metoda Gaussa-Seidla
    std::vector<double> residuumGauss;
    std::vector<double> iterationsGauss;
    auto startGauss = std::chrono::high_resolution_clock::now();
    Matrix rGauss = CalculateGauss(A, b, 1e-9, iterationsGauss, residuumGauss);
    auto stopGauss = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Gaussa-Seidla-"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterationsGauss.size()-1<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopGauss-startGauss).count()<<"ms"<<"\n";
    std::cout<<"Norma bledu rezydualnego: "<<residuumGauss[residuumGauss.size()-1]<<"\n";
    //generowanie wykresu
    plt::figure();
    plt::figure_size(700, 500);
    plt::named_semilogy("Metoda Jacobiego", iterationsJacobi, residuumJacobi);
    plt::named_semilogy("Metoda Gaussa-Seidla", iterationsGauss, residuumGauss);
    plt::title("Porownanie bledu rezydualnego w kolejnych iteracjach (zadanie B)");
    plt::xlabel("Iteracja");
    plt::ylabel("Norma bledu rezydualnego");
    plt::legend();
    plt::save("plots/ZadanieBResiduum.png");
}

void Example2(){
    std::cout<<"\nZadanie C:"<<"\n";
    std::cout<<"N="<<N<<"\n";
    Matrix A(N, N);
    A.CreateDiagonal(3, A2, A3);
    Matrix b(N, 1);
    for(int i=0; i<N; i++){
        b.Set(i, 0, sin(i*(F+1)));
    }
    //metoda Jacobiego
    std::vector<double> residuumJacobi;
    std::vector<double> iterationsJacobi;
    auto startJacobi = std::chrono::high_resolution_clock::now();
    Matrix rJacobi = CalculateJacobi(A, b, 1e-9, iterationsJacobi, residuumJacobi);
    auto stopJacobi = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Jacobiego-"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterationsJacobi.size()-1<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopJacobi-startJacobi).count()<<"ms"<<"\n";
    std::cout<<"Norma bledu rezydualnego: "<<residuumJacobi[residuumJacobi.size()-1]<<"\n";
    //metoda Gaussa-Seidla
    std::vector<double> residuumGauss;
    std::vector<double> iterationsGauss;
    auto startGauss = std::chrono::high_resolution_clock::now();
    Matrix rGauss = CalculateGauss(A, b, 1e-9, iterationsGauss, residuumGauss);
    auto stopGauss = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Gaussa-Seidla-"<<"\n";
    std::cout<<"Liczba iteracji: "<<iterationsGauss.size()-1<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopGauss-startGauss).count()<<"ms"<<"\n";
    std::cout<<"Norma bledu rezydualnego: "<<residuumGauss[residuumGauss.size()-1]<<"\n";
    //generowanie wykresu
    plt::figure();
    plt::figure_size(700, 500);
    plt::named_semilogy("Metoda Jacobiego", iterationsJacobi, residuumJacobi);
    plt::named_semilogy("Metoda Gaussa-Seidla", iterationsGauss, residuumGauss);
    plt::title("Porownanie bledu rezydualnego w kolejnych iteracjach (zadanie C)");
    plt::xlabel("Iteracja");
    plt::ylabel("Norma bledu rezydualnego");
    plt::legend();
    plt::save("plots/ZadanieCResiduum.png");

    std::cout<<"\nZadanie D:"<<"\n";
    std::cout<<"N="<<N<<"\n";
    auto startLU = std::chrono::high_resolution_clock::now();
    Matrix xLU = LUFactorization(A, b);
    auto stopLU = std::chrono::high_resolution_clock::now();
    std::cout<<"-Metoda Faktoryzacji LU-"<<"\n";
    std::cout<<"Czas wykonania: "<<std::chrono::duration_cast<std::chrono::milliseconds>(stopLU-startLU).count()<<"ms"<<"\n";
    std::cout<<"Norma bledu rezydualnego: "<<(A*xLU-b).CalculateNorm()<<"\n";

}

void Example3(){
    std::cout<<"\nZadanie E:"<<"\n";
    std::vector<double> ns = {100, 500, 1000, 2000, 3000};
    std::vector<double> timesJacobi;
    std::vector<double> timesGauss;
    std::vector<double> timesLU;
    std::vector<double> iterationsJacobi;
    std::vector<double> iterationsGauss;
    std::vector<double> rJacobi;
    std::vector<double> rGauss;
    std::vector<double> rLU;

    double epsilon = 1e-9;
    for(double n : ns){
        Matrix A(n, n);
        A.CreateDiagonal(A1, A2, A3);
        Matrix b(n, 1);
        for(int i=0; i<n; i++){
            b.Set(i, 0, sin(i*(F+1)));
        }
        std::vector<double> iJacobi;
        std::vector<double> residuumJacobi;
        auto startJacobi = std::chrono::high_resolution_clock::now();
        Matrix xJacobi = CalculateJacobi(A, b, 1e-9, iJacobi, residuumJacobi);
        auto stopJacobi = std::chrono::high_resolution_clock::now();
        double timeJacobi = std::chrono::duration_cast<std::chrono::milliseconds>(stopJacobi-startJacobi).count();
        timesJacobi.push_back(timeJacobi);
        iterationsJacobi.push_back(iJacobi.size()-1);
        rJacobi.push_back(residuumJacobi[residuumJacobi.size()-1]);
        std::cout<<"Metoda Jacobiego N: "<<n<<" czas: "<<timeJacobi<<"ms iteracje: "<<iJacobi.size()-1<<"\n";

        std::vector<double> iGauss;
        std::vector<double> residuumGauss;
        auto startGauss = std::chrono::high_resolution_clock::now();
        Matrix xGauss = CalculateGauss(A, b, 1e-9, iGauss, residuumGauss);
        auto stopGauss = std::chrono::high_resolution_clock::now();
        double timeGauss = std::chrono::duration_cast<std::chrono::milliseconds>(stopGauss-startGauss).count();
        timesGauss.push_back(timeGauss);
        iterationsGauss.push_back(iGauss.size()-1);
        rGauss.push_back(residuumGauss[residuumGauss.size()-1]);
        std::cout<<"Metoda Gaussa N: "<<n<<" czas: "<<timeGauss<<"ms iteracje: "<<iGauss.size()-1<<"\n";

        auto startLU = std::chrono::high_resolution_clock::now();
        Matrix xLU = LUFactorization(A, b);
        auto stopLU = std::chrono::high_resolution_clock::now();
        double timeLU = std::chrono::duration_cast<std::chrono::milliseconds>(stopLU-startLU).count();
        timesLU.push_back(timeLU);
        rLU.push_back((A*xLU-b).CalculateNorm());
        std::cout<<"Faktoryzacja LU N: "<<n<<" czas: "<<timeLU<<"ms\n";
        std::cout<<"\n";
    }

    plt::figure();
    plt::figure_size(700, 500);
    plt::named_plot("Metoda Jacobiego", ns, timesJacobi);
    plt::named_plot("Metoda Gaussa-Seidla", ns, timesGauss);
    plt::title("Czasy trwania algorytmów iteracyjnych w zależnosci od liczby niewiadomych");
    plt::xlabel("Liczba niewiadomych");
    plt::ylabel("Czas trwania(ms)");
    plt::legend();
    plt::save("plots/ZadanieECzasy.png");

    plt::named_plot("Faktoryzacja LU", ns, timesLU);
    plt::title("Porównanie metod iteracyjnych z bespośrednią faktoryzacją LU");
    plt::legend();
    plt::save("plots/ZadanieECzasy2.png");

    plt::figure();
    plt::named_plot("Metoda Jacobiego",ns, iterationsJacobi);
    plt::named_plot("Metoda Gaussa-Seidla",ns, iterationsGauss); 
    plt::title("Ilość iteracji w zależności od liczby niewiadomych");
    plt::xlabel("Liczba niewiadomych");
    plt::ylabel("Ilość iteracji");
    plt::legend();
    plt::save("plots/ZadanieEIteracje.png");

    plt::figure();
    plt::named_semilogy("Metoda Jacobiego", ns, rJacobi);
    plt::named_semilogy("Metoda Gaussa-Seidla", ns, rGauss);
    plt::named_semilogy("Faktoryzacja LU", ns, rLU);
    plt::title("Norma błędu rezydualnego w zależności od liczby niewiadomych");
    plt::xlabel("Liczba niewiadomych");
    plt::ylabel("Norma residuum");
    plt::legend();
    plt::save("plots/ZadanieEResiduum");
}

int main(){
    //zadanie A i B
    //Example1();
    //zadanie C i D
    //Example2();
    //zadanie E
    Example3();
    
    return 0;
}