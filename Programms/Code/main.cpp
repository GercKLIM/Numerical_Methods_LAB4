#include <iostream>
#include"SolvePDE.h"
int main() {
    std::cout << "Hello, World!" << std::endl;

    // Инициализация первого теста из методички
    double x0_test1 = 0.;
    double X_test1 = 1.;
    double y0_test1 = 0.;
    double Y_test1 = 1.;
    double t0_test1 = 0.;
    double T_test1 = 1.;
    double tau_test1 = 0.001;
    double hx_test1 = 0.2;
    double hy_test1 = 0.2;
    PDEProblem test1(x0_test1, X_test1, y0_test1, Y_test1, t0_test1, T_test1, tau_test1, hx_test1, hy_test1);
    test1.initDeflectionFunc = ([&] (std::vector<double> point) {return 1.;});
    test1.initDeflectionFunc_isSet = true;
    test1.dirichletBoundaryFunc_North = ([&] (std::vector<double> point) {return 1.;});
    test1.dirichletBoundaryFunc_North_isSet = true;
    test1.dirichletBoundaryFunc_South = ([&] (std::vector<double> point) {return 1.;});
    test1.dirichletBoundaryFunc_South_isSet = true;
    test1.dirichletBoundaryFunc_West = ([&] (std::vector<double> point) {return 1.;});
    test1.dirichletBoundaryFunc_West_isSet = true;
    test1.dirichletBoundaryFunc_East = ([&] (std::vector<double> point) {return 1.;});
    test1.dirichletBoundaryFunc_East_isSet = true;
    // Расчёт первого теста
    LongTransScheme(test1, "Test1.txt");
    return 0;
}
