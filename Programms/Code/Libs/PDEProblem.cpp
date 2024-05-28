//
// Created by Иван on 5/24/2024.
//

#include "PDEProblem.h"

PDEProblem::PDEProblem(const double &x0_init, const double &X_init, const double &y0_init, const double &Y_init,
                       const double &t0_init, const double &T_init, const double &tau_init, const double &hx_init,
                       const double &hy_init) {
    x0 = x0_init;
    X = X_init;
    y0 = y0_init;
    Y = Y_init;
    t0 = t0_init;
    T = T_init;
    tau = tau_init;
    hx = hx_init;
    hy = hy_init;
    if(hx == 0 || hy == 0){
        std::cout << "LOG[ERROR] hx or hy is set to 0!!!" << std::endl;
    }
    num_time_steps = static_cast<int>((T-t0) / tau);
    num_x_steps = static_cast<int>((X - x0)/hx);
    num_y_steps = static_cast<int>((Y - y0)/hy);
}
