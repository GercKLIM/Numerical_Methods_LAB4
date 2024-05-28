//
// Created by Иван on 5/24/2024.
//
#include "SolvePDE.h"

//решение СЛАУ методом правой прогонки
//A x_{i-1} - B x_i + C x_{i+1} = -D  (***)
// Если векторы диагоналей исходной системы заданы как A,B,C,D (A - диагональ опд главной, B - главная диагональ, C - диагональ над главной, D- правая часть)
// То для правильного расчёта необходимо передавать A, (-1.)*B, C, (-1.)*D
// Так как прогонка актуальная для системы (***)
std::vector<double> TridiagonalMatrixAlgorithm(
        std::vector<double> a,
        std::vector<double> b,
        std::vector<double> c,
        std::vector<double> d
){
    int n = b.size();
    std::vector<double> alphas({ c[0] / b[0] });
    std::vector<double> betas({d[0] / b[0]});
    for (int i = 1; i < n-1; ++i)
    {
        double denom = b[i] - a[i] * alphas[i - 1];
        alphas.push_back(c[i] / denom);
        betas.push_back((d[i] + a[i] * betas[i-1]) / denom);
    }
    betas.push_back((d[n - 1] + a[n - 1] * betas[n - 2]) / (b[n-1] - a[n-1] * alphas[n - 2]));
    std::vector<double> SolutionX({betas[n-1]});
    for (int i = n - 2; i >= 0; --i)
    {
        SolutionX.push_back(alphas[i] * SolutionX[n - i - 2] + betas[i]);
    }
    reverse(SolutionX.begin(), SolutionX.end());
    return SolutionX;
}

/* Инициализация начального состояния системы (каждой точке пространства ставим в соотвествие величину начального отклонения)
 * args: problem
 * return: initial state
 * */
std::vector<std::vector<double>>  initializeState(const PDEProblem &problem){
    std::vector<std::vector<double>> new_state(problem.num_x_steps+1, std::vector<double>(problem.num_y_steps+1,0));
    double x_i = problem.x0;
    double y_i = problem.y0;

    /*     y0   y1 ...  yM
     *   -------------------
     * x0| u00  u01 ... u0M
     * x1| u10  u11 ... u1M
     *...| ................
     * xN| uN0  uN1 ... uNM
     * ---------------------
     * */

    if(problem.initDeflectionFunc_isSet) {
        for (int i = 0; i <= problem.num_x_steps; ++i) {
            for (int j = 0; j <= problem.num_y_steps; ++j) {
                new_state[i][j] = problem.initDeflectionFunc({x_i, y_i});
                y_i += problem.hy;
            }
            y_i = problem.y0;
            x_i += problem.hx;
        }
    }
    else{
        std::cout << "LOG[WARN]: InitDeflectionFunc was not set!!!" << std::endl;
    }

    return new_state;
}

// Функция F из методички
//u__m -- u -- u__p --->y
double F(const double& u, const double& u__m, const double& u__p, const double& x, const double& y, const PDEProblem& problem){
    double res = 2/problem.tau * u + (u__p - 2*u + u__m)/(problem.hy * problem.hy);
    if(problem.extForcesFunction_isSet){
        res += problem.extForcesFunction({x,y});
    }
    return res;
}

// Функция F с шапочкой (из фольги хи-хи) из методички
//u_m_ -- u -- u_p_ --->x
double F_hat(const double& u, const double& u_m_, const double& u_p_, const double& x, const double& y, const PDEProblem& problem){
    double res = 2/problem.tau * u + (u_p_ - 2*u + u_m_)/(problem.hx * problem.hx);
    if(problem.extForcesFunction_isSet){
        res += problem.extForcesFunction({x,y});
    }
    return res;
}

bool LongTransScheme(const PDEProblem &problem, const string &filename) {

    // Инициализация начального состояния
    std::vector<std::vector<double>> state_k = initializeState(problem);
    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open()) {
        double t_i = problem.t0;
        double x_i = problem.x0;
        double y_i = problem.y0;
        double tau = problem.tau;
        double half_tau = tau/2;
        double hx = problem.hx;
        double hy = problem.hy;
        std::vector<std::vector<double>> state_kp(state_k);
        fpoints << t_i << endl;
        write2DVectorToFile(fpoints, state_k);

        // Utility vectors
        std::vector<double> Axs(problem.num_x_steps+1, 0);
        std::vector<double> Bxs(problem.num_x_steps+1, 0);
        std::vector<double> Cxs(problem.num_x_steps+1, 0);
        std::vector<double> Dxs(problem.num_x_steps+1, 0);
        std::vector<double> Ays(problem.num_x_steps+1, 0);
        std::vector<double> Bys(problem.num_y_steps+1, 0);
        std::vector<double> Cys(problem.num_y_steps+1, 0);
        std::vector<double> Dys(problem.num_y_steps+1, 0);

        for(int j = 0; j < problem.num_time_steps; ++j){
            // слой ,,половинный'' первый (k+1/2)
            t_i += half_tau;

            // Calculation across X-axis
            for(int i2 = 1; i2 < problem.num_y_steps; ++i2) {
                y_i += hy;

                // Г.У. Запад
                if (problem.dirichletBoundaryFunc_West_isSet) {
                    Axs[0] = 0.;
                    Bxs[0] = 1.;
                    Cxs[0] = 0.;
                    Dxs[0] = problem.dirichletBoundaryFunc_West({problem.x0, y_i});
                } else if (problem.neymanBoundaryFunc_West_isSet) {
                    /* аппроксимация второго рода*/
                }

                // Г.У. Восток
                if (problem.dirichletBoundaryFunc_East_isSet) {
                    Axs[problem.num_x_steps] = 0.;
                    Bxs[problem.num_x_steps] = 1.;
                    Cxs[problem.num_x_steps] = 0.;
                    Dxs[problem.num_x_steps] = problem.dirichletBoundaryFunc_East({problem.X, y_i});
                } else if (problem.neymanBoundaryFunc_East_isSet) {
                    /* аппроксимация второго рода*/
                }

                for(int i1 = 1; i1 < problem.num_x_steps; ++i1){
                    x_i += hx;
                    Axs[i1] = 1/(hx*hx);
                    Bxs[i1] = -2*( 1/(hx*hx) + 1/tau );
                    Cxs[i1] = 1/(hx*hx);
                    Dxs[i1] = -F(state_k[i1][i2], state_k[i1][i2-1],state_k[i1][i2+1], x_i, y_i, problem);
                }

                std::vector<double> x_line_sol = TridiagonalMatrixAlgorithm(Axs, Bxs, Cxs, Dxs);
                for(int i1 = 0; i1<=problem.num_x_steps; ++i1){
                    state_kp[i1][i2] = x_line_sol[i1];
                }
            }
            x_i = problem.x0;
            y_i = problem.y0;
            // слой ,,половинный'' второй (k+1)
            t_i += half_tau;

            // Calculation across Y-axis
            for(int i1 = 1; i1 < problem.num_x_steps; ++i1) {
                x_i += hx;

                // Г.У. Север
                if (problem.dirichletBoundaryFunc_North_isSet) {
                    Ays[problem.num_y_steps] = 0.;
                    Bys[problem.num_y_steps] = 1.;
                    Cys[problem.num_y_steps] = 0.;
                    Dys[problem.num_y_steps] = problem.dirichletBoundaryFunc_North({x_i, problem.Y});
                } else if (problem.neymanBoundaryFunc_North_isSet) {
                    /* аппроксимация второго рода*/
                }
                // Г.У. Юг
                if (problem.dirichletBoundaryFunc_South_isSet) {
                    Ays[0] = 0.;
                    Bys[0] = 1.;
                    Cys[0] = 0.;
                    Dys[0] = problem.dirichletBoundaryFunc_South({x_i, problem.y0});
                } else if (problem.neymanBoundaryFunc_South_isSet) {
                    /* аппроксимация второго рода*/
                }

                for(int i2 = 1; i2 < problem.num_y_steps; ++i2){
                    Ays[i2] = 1/(hy*hy);
                    Bys[i2] = -2*( 1/(hy*hy) + 1/tau );
                    Cys[i2] = 1/(hy*hy);
                    Dys[i2] = -F_hat(state_kp[i1][i2], state_kp[i1-1][i2],state_kp[i1+1][i2], x_i, y_i, problem);
                }

                std::vector<double> y_line_sol = TridiagonalMatrixAlgorithm(Ays, Bys, Cys, Dys);
                for(int i2 = 0; i2<=problem.num_y_steps; ++i2){
                    state_kp[i1][i2] = y_line_sol[i2];
                }
            }
            x_i = problem.x0;
            y_i = problem.y0;
            fpoints << t_i << endl;
            write2DVectorToFile(fpoints, state_kp);
        }
        return true;
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
/*
    // Инициализация начального состояния
    //std::vector<double> state_0 = init_state(num_space_steps, u_0); //TODO: расширить init_state
    std::vector<double> state_0 = init_state(num_space_steps+1, h, test);
    std::vector<double> As(num_space_steps+1, 0);
    std::vector<double> Cs(num_space_steps+1, 0);
    std::vector<double> Bs(num_space_steps+1, 0);
    std::vector<double> Fs(num_space_steps+1, 0);

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        double t_i = t_0;
        std::vector<double> state_i = state_0;
        writeVectorToFile(fpoints, t_i, state_i);
        double x_i = x_0;

        // Эволюция системы во времени
        for(int j = 0; j <= num_time_steps; ++j) {
            t_i += tau;

            // Граничные условия слева

            // 1-го рода
            if(!test.G_left_type){
                Cs[0] = 1.;
                Bs[0] = 0.;
                As[0] = 0.;
                Fs[0] = test.G_left(t_i);
            }

                // 2-го рода
            else {
                double a1 = a(test.K_ptr, x_0+h, x_0);
                double w1 = w(a1, state_i[1], state_i[0], h);
                double kappa = sigma*a1/h / (c*rho*h/(2*tau)+sigma*a1/h);
                double mu = (c*rho*state_i[0]*h/(2*tau)+sigma*test.G_left(t_i)+(1-sigma)*(test.G_left(t_i-tau)+w1))/(c*rho*h/(2*tau)+sigma*a1/h);
                Cs[0] = 1.;
                Bs[0] = kappa;
                As[0] = 0;
                Fs[0] = mu;

                if(max_a < fabs(a1))
                {
                    max_a = fabs(a1);
                }
            }

            // Граничные условия справа
            // 1-го рода
            if(!test.G_right_type){
                Bs[num_space_steps] = 0.;
                As[num_space_steps] = 0.;
                Cs[num_space_steps] = 1.;
                Fs[num_space_steps] = test.G_right(t_i);
            }

                // 2-го рода
            else{
                double am = a(test.K_ptr, X, X-h);
                double wn = w(am, state_i[num_space_steps], state_i[num_space_steps-1], h);
                double denom = c * rho * h / (2 * tau) + sigma * am / h;
                double kappa = sigma * am /h / denom;
                double mu = (c * rho * state_i[num_space_steps] * h / (2 * tau) + sigma * test.G_right(t_i) + (1 - sigma) * (test.G_right(t_i-tau) - wn)) / denom;
                Cs[num_space_steps] = 1.;
                Bs[num_space_steps] = 0.;
                As[num_space_steps] = kappa;
                Fs[num_space_steps] = mu;
                if(max_a < fabs(am))
                {
                    max_a = fabs(am);
                }
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                double a_i = a(test.K_ptr, x_i, x_i - h);
                double a_ip = a(test.K_ptr, x_i + h, x_i);
                if(max_a < fabs(a_i))
                {
                    max_a = fabs(a_i);
                }
                if(max_a < fabs(a_ip))
                {
                    max_a = fabs(a_ip);
                }
                As[i] = sigma / h * a_i;
                Bs[i] = sigma / h * a_ip;
                Cs[i] = As[i] + Bs[i] + c * rho * h / tau;
                Fs[i] = c * rho * h / tau * state_i[i] +
                        (1 - sigma) * (w(a_ip, state_i[i + 1], state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h));
            }

            // Получение нового состояния системы
            // A - C + B = - F (не домножаем векторы на -1, так как уже считали домноженные)
            state_i = TridiagonalMatrixAlgorithm(As, Cs, Bs, Fs);

            // Запись в файл
            writeVectorToFile(fpoints, t_i, state_i);
        }
        fpoints.close();
        std::cout << "!!!!!MAX A = " << max_a << std::endl;
        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
    */
};
