//
// Created by Иван on 5/24/2024.
//

#ifndef CODE_PDEPROBLEM_H
#define CODE_PDEPROBLEM_H

#include <cmath>
#include <iostream>
#include <functional>
#include "Libs/algebra.h"

class PDEProblem {
public:
    double x0;  // начало отсчёта по x
    double X;   // координата правого конца по x
    double y0;  // начало отсчёта по y
    double Y;   // кордината правого конца по Y
    double t0;  // время отсчёта
    double T;   // время окончания
    double tau; // шаг по времени
    double hx;  // шаг по x-координате
    double hy;  // шаг по y-координате
    int num_time_steps; // количество шагов по времени
    int num_x_steps; // количество шагов по пространству x
    int num_y_steps; // количество шагов по пространству y

    PDEProblem(const double &x0_init, const double &X_init, const double &y0_init, const double &Y_init,
               const double &t0_init, const double &T_init, const double &tau_init, const double &hx_init,
               const double &hy_init);

    /*
 *         North
 *           ^
 *       ____|____
 *      |    |    |
 * West |- - * - -|->East
 *      |   0|    |
 *      |____|____|
 *         South
 * */

    // Отклонение точки в нулевой момент времени (U_0(x, y, 0))
    std::function<double(std::vector<double>)> initDeflectionFunc;
    bool  initDeflectionFunc_isSet = false;

    // Начальное условие Дирихле для северной границы \xi(x,y,t)
    std::function<double(std::vector<double>)> dirichletBoundaryFunc_North;
    bool  dirichletBoundaryFunc_North_isSet = false;

    // Начальное условие Дирихле для южной границы \xi(x,y,t)
    std::function<double(std::vector<double>)> dirichletBoundaryFunc_South;
    bool  dirichletBoundaryFunc_South_isSet = false;

    // Начальное условие Дирихле для западной границы \xi(x,y,t)
    std::function<double(std::vector<double>)> dirichletBoundaryFunc_West;
    bool  dirichletBoundaryFunc_West_isSet = false;

    // Начальное условие Дирихле для восточной границы \xi(x,y,t)
    std::function<double(std::vector<double>)> dirichletBoundaryFunc_East;
    bool  dirichletBoundaryFunc_East_isSet = false;

    // Начальное условие Неймана для северной границы \psi(x,y,t)
    std::function<double(std::vector<double>)> neymanBoundaryFunc_North;
    bool  neymanBoundaryFunc_North_isSet = false;

    // Начальное условие Неймана для южной границы \psi(x,y,t)
    std::function<double(std::vector<double>)> neymanBoundaryFunc_South;
    bool  neymanBoundaryFunc_South_isSet = false;

    // Начальное условие Неймана для западной границы \psi(x,y,t)
    std::function<double(std::vector<double>)> neymanBoundaryFunc_West;
    bool  neymanBoundaryFunc_West_isSet = false;

    // Начальное условие Неймана для восточной границы \psi(x,y,t)
    std::function<double(std::vector<double>)> neymanBoundaryFunc_East;
    bool  neymanBoundaryFunc_East_isSet = false;

    // Воздействие внешних сил (F(x,y)) (если == 0, то уравнение - однородное)
    std::function<double(std::vector<double>)> extForcesFunction;
};

#endif //CODE_PDEPROBLEM_H


