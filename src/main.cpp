#include <iostream>
#include <vector>
#include <matplotlibcpp.h>
#include "matrix.h"

#define PLOT_DIR "plots/"

namespace plt = matplotlibcpp;

int main(){
    Matrix* m = new Matrix(2, 2, {{2.0, 3.0},{8.0, 76.0}});
    Matrix* d = new Matrix(2, 2, {{1.0, 1.0},{1.0, 1.0}});
    std::cout<<m->ToString()<<"\n";
    Matrix* c = m->Add(*d);
    std::cout<<m->Add(*d)->ToString()<<"\n";
    return 0;
}