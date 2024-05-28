//
// Created by Иван on 5/24/2024.
//
#include "SolvePDE.h"
/* Инициализация начального состояния системы (каждой точке пространства ставим в соотвествие величину начального отклонения)
 * args: problem
 * return: initial state
 * */
std::vector<std::vector<double>>  initializeState(const PDEProblem &problem){
    std::vector<std::vector<double>> new_state(problem.num_x_steps+1, std::vector<double>(problem.num_y_steps+1,0));
    double x_i = problem.x0;
    double y_i = problem.y0;
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



bool LongTransScheme(const PDEProblem &problem, const string &filename) {

    // Инициализация начального состояния
    std::vector<std::vector<double>> state_j = initializeState(problem);
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
        std::vector<std::vector<double>> state_jp(state_j);
        fpoints << t_i << endl;
        write2DVectorToFile(fpoints, state_j);
        for(int j = 0; j < problem.num_time_steps; ++j){
            t_i += half_tau;

            // Calculation across X-axis

            t_i += half_tau;

            // Calculation across Y-axis

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
