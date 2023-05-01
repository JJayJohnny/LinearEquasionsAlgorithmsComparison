#include <iostream>
#include <vector>
#include <matplotlibcpp.h>
#include "matrix.h"

#define PLOT_DIR "plots/"

namespace plt = matplotlibcpp;

int main(){
    Matrix m(5, 2, {{5.0, 1.0},{-3.0, 3.0},{4.0, -1.0},{-1.0, 2.0},{0.0, 7.0}});
    Matrix d(2, 7, {{1.0, 7.0, 0.0, -1.0, 9.0, 0.0, 5.0},{2.0, -6.0, -1.0, 10.0, 2.0, -3.0, 1.0}});
    std::cout<<m.ToString()<<"\n";
    std::cout<<(m+1).ToString()<<"\n";
    return 0;
}