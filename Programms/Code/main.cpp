#include <iostream>
#include"SolvePDE.h"


void make_data_for_tables_test1(){

    // Инициализация первого теста из методички
    double x0_test1 = 0.;
    double X_test1 = 1.;
    double y0_test1 = 0.;
    double Y_test1 = 1.;
    double t0_test1 = 0.;
    double T_test1 = 10.;
    double tau = 10e-4;

    double h0_x = 0.2;
    double h0_y = 0.2;

    for (int i = 1; i < 6; i++ ){

        // Теструем на тесте 1 из методички
        PDEProblem test1(x0_test1, X_test1, y0_test1, Y_test1, t0_test1, T_test1, tau,  h0_x, h0_y);
        test1.initDeflectionFunc = ([&] (std::vector<double> point) {return 10.;});
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
        tau = tau / 2.;

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
    double T_test2 = 10.;
    double tau = 10e-8;
    double h0_x = 0.2;
    double h0_y = 0.2;


    for (int i = 1; i < 6; i++ ){

        // Теструем на тесте 2 из методички
        PDEProblem test2(x0_test2, X_test2, y0_test2, Y_test2, t0_test2, T_test2, tau, h0_x, h0_y);
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
        tau = tau / 2.;
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
    double hx_test1 = 0.2;
    double hy_test1 = 0.2;
    double tau_test1 = std::sqrt(hx_test1*hx_test1 + hy_test1*hy_test1)/std::sqrt(1./(X_test1*X_test1) + 1./(Y_test1*Y_test1))/M_PI;
    PDEProblem test1(x0_test1, X_test1, y0_test1, Y_test1, t0_test1, T_test1, tau_test1, hx_test1, hy_test1);
    test1.initDeflectionFunc = ([&] (std::vector<double> point) {return 0.9;});
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
}

void test2(){
    // Инициализация второго теста из методички
    double x0_test2 = 0.;
    double X_test2 = 1.;
    double y0_test2 = 0.;
    double Y_test2 = 1.;
    double t0_test2 = 0.;
    double T_test2 = 10.;
    //double tau_test2 = 0.00001;
    double hx_test2 = 0.2;
    double hy_test2 = 0.2;
    double tau_test2 = std::sqrt(hx_test2*hx_test2 + hy_test2*hy_test2)/std::sqrt(1./(X_test2*X_test2) + 1./(Y_test2*Y_test2))/M_PI;

    PDEProblem test2(x0_test2, X_test2, y0_test2, Y_test2, t0_test2, T_test2, tau_test2, hx_test2, hy_test2);
    test2.initDeflectionFunc = ([&] (std::vector<double> point) {return -1.;});
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


void test3(){
    // Инициализация второго теста из методички
    double x0_test3 = 0.;
    double X_test3 = 1.;
    double y0_test3 = 0.;
    double Y_test3 = 1.;
    double t0_test3 = 0.;
    double T_test3 = 3.;
    double hx_test3 = 0.1;
    double hy_test3 = 0.1;
    double tau_test3 = std::sqrt(hx_test3*hx_test3 + hy_test3*hy_test3)/std::sqrt(1./(X_test3*X_test3) + 1./(Y_test3*Y_test3))/M_PI;

    PDEProblem test3(x0_test3, X_test3, y0_test3, Y_test3, t0_test3, T_test3, tau_test3, hx_test3, hy_test3);
    test3.initDeflectionFunc = ([&] (std::vector<double> point) {return 1. + point[1];});
    test3.initDeflectionFunc_isSet = true;
    test3.extForcesFunction = ([&] (std::vector<double> point) {return 4.;});
    test3.extForcesFunction_isSet = true;
    test3.neymanBoundaryFunc_West = ([&] (std::vector<double> point) {return 0.;});
    test3.neymanBoundaryFunc_West_isSet = true;
    test3.neymanBoundaryFunc_East = ([&] (std::vector<double> point) {return 2.;});
    test3.neymanBoundaryFunc_East_isSet = true;
    test3.dirichletBoundaryFunc_South = ([&] (std::vector<double> point) {return point[0] * point[0];});
    test3.dirichletBoundaryFunc_South_isSet = true;
    test3.dirichletBoundaryFunc_North = ([&] (std::vector<double> point) {return 1. + point[0] * point[0];});
    test3.dirichletBoundaryFunc_North_isSet = true;
    // Расчёт третьего тестаё
    LongTransScheme(test3, "Test3.txt");
}

void test4(){
    // Инициализация четвёртого теста
    double x0_test4 = 0.;
    double X_test4 = 1.;
    double y0_test4 = 0.;
    double Y_test4 = 1.;
    double t0_test4 = 0.;
    double T_test4 = 10.;
    //double tau_test2 = 0.00001;
    double hx_test4 = 0.2;
    double hy_test4 = 0.2;
    double tau_test4 = std::sqrt(hx_test4*hx_test4 + hy_test4*hy_test4)/std::sqrt(1./(X_test4*X_test4) + 1./(Y_test4*Y_test4))/M_PI;

    PDEProblem test4(x0_test4, X_test4, y0_test4, Y_test4, t0_test4, T_test4, tau_test4, hx_test4, hy_test4);
    test4.initDeflectionFunc = ([&] (std::vector<double> point) {return -1.;});
    test4.initDeflectionFunc_isSet = true;
    test4.dirichletBoundaryFunc_North = ([&] (std::vector<double> point) {return 1.+point[0];});
    test4.dirichletBoundaryFunc_North_isSet = true;
    test4.dirichletBoundaryFunc_South = ([&] (std::vector<double> point) {return 1.+point[0];});
    test4.dirichletBoundaryFunc_South_isSet = true;
    test4.neymanBoundaryFunc_West = ([&] (std::vector<double> point) {return -1.;});
    test4.neymanBoundaryFunc_West_isSet = true;
    test4.neymanBoundaryFunc_East = ([&] (std::vector<double> point) {return 1.;});
    test4.neymanBoundaryFunc_East_isSet = true;
    // Расчёт второго теста
    LongTransScheme(test4, "Test4.txt");
}


int main() {
    //make_data_for_tables_test1();
    //make_data_for_tables_test2();
    //test1();
     //test2();
    //test3();
    test4();
    return 0;
}
