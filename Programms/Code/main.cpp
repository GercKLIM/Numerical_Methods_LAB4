#include <iostream>
#include"SolvePDE.h"


void make_data_for_tables_test1(){

    // Инициализация первого теста из методички
    double x0_test1 = 0.;
    double X_test1 = 1.;
    double y0_test1 = 0.;
    double Y_test1 = 1.;
    double t0_test1 = 0.;
    double T_test1 = 1.;
    double tau_test1 = 0.001;

    double h0_x = 0.2;
    double h0_y = 0.2;

    for (int i = 1; i < 6; i++ ){

        // Теструем на тесте 1 из методички
        PDEProblem test1(x0_test1, X_test1, y0_test1, Y_test1, t0_test1, T_test1, tau_test1,  h0_x, h0_y);
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

        // Изменение шага
        h0_x = h0_x / 2.;
        h0_y = h0_y / 2.;

        // Расчёт первого теста
        LongTransScheme_for_tables(test1, "data_for_tables/test1/test1_" + to_string(i) + ".txt");
    }
}


void make_data_for_tables_test2(){

    // Инициализация второго теста из методички
    double x0_test2 = 0.;
    double X_test2 = 1.;
    double y0_test2 = 0.;
    double Y_test2 = 1.;
    double t0_test2 = 0.;
    double T_test2 = 1.;
    double tau_test2 = 0.001;
    double h0_x = 0.2;
    double h0_y = 0.2;


    for (int i = 1; i < 6; i++ ){

        // Теструем на тесте 2 из методички
        PDEProblem test2(x0_test2, X_test2, y0_test2, Y_test2, t0_test2, T_test2, tau_test2, h0_x, h0_y);
        test2.initDeflectionFunc = ([&] (std::vector<double> point) {return 1. + point[1];});
        test2.initDeflectionFunc_isSet = true;
        test2.neymanBoundaryFunc_North = ([&] (std::vector<double> point) {return 1.;});
        test2.neymanBoundaryFunc_North_isSet = true;
        test2.neymanBoundaryFunc_South = ([&] (std::vector<double> point) {return -1.;});
        test2.neymanBoundaryFunc_South_isSet = true;
        test2.dirichletBoundaryFunc_West = ([&] (std::vector<double> point) {return 1. + point[1];});
        test2.dirichletBoundaryFunc_West_isSet = true;
        test2.dirichletBoundaryFunc_East = ([&] (std::vector<double> point) {return 1. + point[1];});
        test2.dirichletBoundaryFunc_East_isSet = true;

        // Изменение шага
        h0_x = h0_x / 2.;
        h0_y = h0_y / 2.;

        // Расчёт первого теста
        LongTransScheme_for_tables(test2, "data_for_tables/test2/test2_" + to_string(i) + ".txt");
    }


}


void test1(){

    // Инициализация первого теста из методички
    double x0_test1 = 0.;
    double X_test1 = 1.;
    double y0_test1 = 0.;
    double Y_test1 = 1.;
    double t0_test1 = 0.;
    double T_test1 = 1.;
    double tau_test1 = 0.1;
    double hx_test1 = 0.02;
    double hy_test1 = 0.02;
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
    LongTransScheme_for_tables(test1, "Test2.txt");
}

void test2(){
    // Инициализация второго теста из методички
    double x0_test2 = 0.;
    double X_test2 = 1.;
    double y0_test2 = 0.;
    double Y_test2 = 1.;
    double t0_test2 = 0.;
    double T_test2 = 10.;
    double tau_test2 = 0.00001;
    double hx_test2 = 0.1;
    double hy_test2 = 0.1;
    PDEProblem test2(x0_test2, X_test2, y0_test2, Y_test2, t0_test2, T_test2, tau_test2, hx_test2, hy_test2);
    test2.initDeflectionFunc = ([&] (std::vector<double> point) {return 1. + point[1];});
    test2.initDeflectionFunc_isSet = true;
    test2.neymanBoundaryFunc_North = ([&] (std::vector<double> point) {return 1.;});
    test2.neymanBoundaryFunc_North_isSet = true;
    test2.neymanBoundaryFunc_South = ([&] (std::vector<double> point) {return -1.;});
    test2.neymanBoundaryFunc_South_isSet = true;
    test2.dirichletBoundaryFunc_West = ([&] (std::vector<double> point) {return 1. + point[1];});
    test2.dirichletBoundaryFunc_West_isSet = true;
    test2.dirichletBoundaryFunc_East = ([&] (std::vector<double> point) {return 1. + point[1];});
    test2.dirichletBoundaryFunc_East_isSet = true;
    // Расчёт второго теста
    LongTransScheme(test2, "Test2.txt");
}



int main() {
    //make_data_for_tables_test1();
    make_data_for_tables_test2();

    return 0;
}
