//
// Created by Иван on 5/24/2024.
//
#include "SolvePDE.h"
//
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

    /*                   WEST
     *               y0   y1 ...  yM
     *          -------------------
     *          x0| u00  u01 ... u0M
     *          x1| u10  u11 ... u1M    NORTH
     * SOUTH   ...| ................
     *          xN| uN0  uN1 ... uNM
     *          ---------------------
     *                  EAST
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

    if(problem.dirichletBoundaryFunc_North_isSet){
        new_state[0][problem.num_y_steps] = problem.dirichletBoundaryFunc_North({problem.x0, problem.Y});
        new_state[problem.num_x_steps][problem.num_y_steps] = problem.dirichletBoundaryFunc_North({problem.X, problem.Y});
    }
    if(problem.dirichletBoundaryFunc_South_isSet){
        new_state[0][0] = problem.dirichletBoundaryFunc_South({problem.x0, problem.y0});
        new_state[problem.num_x_steps][0] = problem.dirichletBoundaryFunc_South({problem.X, problem.y0});
    }
    if(problem.dirichletBoundaryFunc_West_isSet){
        new_state[0][0] = problem.dirichletBoundaryFunc_West({problem.x0, problem.y0});
        new_state[0][problem.num_y_steps] = problem.dirichletBoundaryFunc_West({problem.x0, problem.Y});
    }
    if(problem.dirichletBoundaryFunc_East_isSet){
        new_state[problem.num_x_steps][0] = problem.dirichletBoundaryFunc_East({problem.X, problem.y0});
        new_state[problem.num_x_steps][problem.num_y_steps] = problem.dirichletBoundaryFunc_East({problem.X, problem.Y});
    }
    return new_state;
}

// Функция F из методички
//u__m -- u -- u__p --->y
double F(const double& u, const double& u__m, const double& u__p, const double& x, const double& y, const PDEProblem& problem){
    double res = 2/problem.tau * u + (u__p - 2*u + u__m)/(problem.hy * problem.hy);
    if(problem.extForcesFunction_isSet){
        res -= problem.extForcesFunction({x,y});
    }
    return res;
}

// Функция F с шапочкой (из фольги хи-хи) из методички
//u_m_ -- u -- u_p_ --->x
double F_hat(const double& u, const double& u_m_, const double& u_p_, const double& x, const double& y, const PDEProblem& problem){
    double res = 2/problem.tau * u + (u_p_ - 2*u + u_m_)/(problem.hx * problem.hx);
    if(problem.extForcesFunction_isSet){
        res -= problem.extForcesFunction({x,y});
    }
    return res;
}

bool Error_calc(const double& eps, const std::vector<std::vector<double>>& state1, const std::vector<std::vector<double>>& state2, const std::vector<std::vector<double>>& state3){

    double cur_err = 0.;
    double n1 = 0.;
    double n2 = 0.;

    for(int i = 1; i < state1.size()-1; ++i){
        for(int j = 1; j < state1[0].size()-1; ++j){
            cur_err = std::fabs(state3[i][j] - state2[i][j]);
            if(cur_err > n1){
                n1 = cur_err;
            }
            cur_err = std::fabs(state2[i][j] - state1[i][j]);
            if(cur_err > n2){
                n2 = cur_err;
            }
        }
    }
    double nu = n1/n2;
    return n1 >= eps*(1-nu);
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
        double hx = problem.hx;
        double hy = problem.hy;
        double tau = problem.tau;
        //double tau = std::sqrt(hx*hx + hy*hy)/std::sqrt(1./(problem.X*problem.X) + 1./(problem.Y*problem.Y))/M_PI;
        std::cout<<tau<<std::endl;
        double half_tau = tau/2;

        std::vector<std::vector<double>> state_kp(state_k);
        std::vector<std::vector<double>> state_kph(state_k);
        fpoints << t_i << endl;
        write2DVectorToFile(fpoints, state_k);

        // Utility vectors
        std::vector<double> Axs(problem.num_x_steps+1, 0.);
        std::vector<double> Bxs(problem.num_x_steps+1, 0.);
        std::vector<double> Cxs(problem.num_x_steps+1, 0.);
        std::vector<double> Dxs(problem.num_x_steps+1, 0.);
        std::vector<double> Ays(problem.num_x_steps+1, 0.);
        std::vector<double> Bys(problem.num_y_steps+1, 0.);
        std::vector<double> Cys(problem.num_y_steps+1, 0.);
        std::vector<double> Dys(problem.num_y_steps+1, 0.);

        // Задаём внешние силы (если нету, то ноль будем прибавлять в коэфах)
        std::function<double(std::vector<double>)> f = ([&] (const std::vector<double>& point) {return  0.;});
        if(problem.extForcesFunction_isSet){
            f = ([&] (const std::vector<double>& point) {return  -problem.extForcesFunction(point);});
        }

        // Максмальное число итераций
        int max_allowed_its = 5000;
        int cur_its = 0;
        do{
            // слой ,,половинный'' первый (k+1/2)
            t_i += half_tau;
            //std::cout << "t_i = " << t_i  << "; cur_its = " << cur_its <<std::endl;
            x_i = problem.x0;
            y_i = problem.y0;

            state_k = state_kp;
            state_kph = state_k;

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
                    Axs[0] = 0.;
                    Bxs[0] = -1./*-(1./tau + 1./(hx*hx))*/;
                    Cxs[0] = -1. / (hx * hx / tau + 1.)/*-1./(hx*hx)*/;
                    Dxs[0] = -(state_k[0][i2] * hx / tau + problem.neymanBoundaryFunc_West({x_i,y_i})) / (hx / tau + 1. / hx)/*-( (1./(2.*hy*hy)) * (state_k[0][i2+1] - 2.*state_k[0][i2] + state_k[0][i2-1]) + (1./tau)*state_k[0][i2] + (1./hx) * problem.neymanBoundaryFunc_West({problem.x0,y_i}) + (0.5)*f({problem.x0,y_i}))*/;
                }

                // Г.У. Восток
                if (problem.dirichletBoundaryFunc_East_isSet) {
                    Axs[problem.num_x_steps] = 0.;
                    Bxs[problem.num_x_steps] = 1.;
                    Cxs[problem.num_x_steps] = 0.;
                    Dxs[problem.num_x_steps] = problem.dirichletBoundaryFunc_East({problem.X, y_i});
                } else if (problem.neymanBoundaryFunc_East_isSet) {
                    /* аппроксимация второго рода*/
                    Axs[problem.num_x_steps] = -1. / (hx * hx / tau + 1.)/*-1./(hx*hx)*/;
                    Bxs[problem.num_x_steps] = -1./*-(1./tau + 1./(hx*hx))*/;
                    Cxs[problem.num_x_steps] = 0.;
                    Dxs[problem.num_x_steps] = -(state_k[problem.num_x_steps][i2] * hx / tau + problem.neymanBoundaryFunc_East({x_i,y_i})) / (hx / tau + 1. / hx) /*-((1./(2*hy*hy))*(state_k[problem.num_x_steps][i2+1] - 2.*state_k[problem.num_x_steps][i2] + state_k[problem.num_x_steps][i2-1]) + (1./hx)*problem.neymanBoundaryFunc_East({problem.X, y_i}) + 1/tau*state_k[problem.num_x_steps][i2] + (1./4.)*(f({problem.X,y_i})+f({problem.X-hx/2.,y_i})))*/;
                }

                for(int i1 = 1; i1 < problem.num_x_steps; ++i1){
                    x_i += hx;
                    Axs[i1] = 1/(hx*hx);
                    Bxs[i1] = 2*( 1/(hx*hx) + 1/tau );
                    Cxs[i1] = 1/(hx*hx);
                    Dxs[i1] = F(state_k[i1][i2], state_k[i1][i2-1],state_k[i1][i2+1], x_i, y_i, problem);
                }

                std::vector<double> x_line_sol = TridiagonalMatrixAlgorithm(Axs, Bxs, Cxs, Dxs);
                for(int i1 = 0; i1<=problem.num_x_steps; ++i1){
                    state_kph[i1][i2] = x_line_sol[i1];
                    //out(state_k[i1]);
                }
                state_kp[0][i2] = state_kph[0][i2];
                state_kp[problem.num_x_steps][i2] = state_kph[problem.num_x_steps][i2];
                //std::cout<<"tridig sol: "<<std::endl;
                //out(x_line_sol);
                //std::cout<<"---"<<std::endl;
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
                    Ays[problem.num_y_steps] = -1. / (hy * hy / tau + 1)/*-1./(hy*hy)*/;
                    Bys[problem.num_y_steps] = -1./*-(1./tau + 1./(hy*hy))*/;
                    Cys[problem.num_y_steps] = 0.;
                    Dys[problem.num_y_steps] = -(state_kph[i1][problem.num_y_steps] * hy / tau + problem.neymanBoundaryFunc_North({x_i, y_i})) / (hy / tau + 1. / hy)/*-((1./(2.*hx*hx))*(state_kph[i1+1][problem.num_y_steps] - 2.*state_kph[i1][problem.num_y_steps] + state_kph[i1-1][problem.num_y_steps]) + (1./hy)*problem.neymanBoundaryFunc_North({x_i, problem.Y}) + 1./tau*state_kph[i1][problem.num_y_steps] + (0.5)*f({x_i,problem.Y}))*/;
                }
                // Г.У. Юг
                if (problem.dirichletBoundaryFunc_South_isSet) {
                    Ays[0] = 0.;
                    Bys[0] = 1.;
                    Cys[0] = 0.;
                    Dys[0] = problem.dirichletBoundaryFunc_South({x_i, problem.y0});
                } else if (problem.neymanBoundaryFunc_South_isSet) {
                    /* аппроксимация второго рода*/
                    Ays[0] = 0.;
                    Bys[0] = -1./*-(1./tau + 1./(hy*hy))*/;
                    Cys[0] = -1. / (hy * hy / tau + 1.)/*-1./(hy*hy)*/;
                    Dys[0] = -(state_kph[i1][0] * hy / tau + problem.neymanBoundaryFunc_South({x_i, y_i})) / (hy / tau + 1. / hy)/*-( (1./(2.*hx*hx)) * (state_kp[i1+1][0] - 2*state_kp[i1][0] + state_kph[i1-1][0]) + (1./tau)*state_kp[i1][0] + 1./tau*state_kph[i1][0] + (1./hy) * problem.neymanBoundaryFunc_South({x_i,problem.y0}) + (0.5)*f({x_i,problem.y0}))*/;
                }

                for(int i2 = 1; i2 < problem.num_y_steps; ++i2){
                    Ays[i2] = 1/(hy*hy);
                    Bys[i2] = 2*( 1/(hy*hy) + 1/tau );
                    Cys[i2] = 1/(hy*hy);
                    Dys[i2] = F_hat(state_kph[i1][i2], state_kph[i1-1][i2],state_kph[i1+1][i2], x_i, y_i, problem);
                }

                std::vector<double> y_line_sol = TridiagonalMatrixAlgorithm(Ays, Bys, Cys, Dys);
                for(int i2 = 0; i2<=problem.num_y_steps; ++i2) {
                    state_kp[i1][i2] = y_line_sol[i2];
                }

                //std::cout<<"tridig sol: "<<std::endl;
                //out(y_line_sol);
                //std::cout<<"---"<<std::endl;
            }
            x_i = problem.x0;
            y_i = problem.y0;
//            fpoints << t_i << endl;
//           write2DVectorToFile(fpoints, state_kp);
            //state_k.swap(state_kp);
            ++cur_its;
        }while(Error_calc(1e-17, state_k, state_kph, state_kp) && cur_its < max_allowed_its/*norm1(state_k + (-1)*state_kp) > 1e-17 && cur_its < max_allowed_its*/ );
        std::cout << cur_its << std::endl;
        fpoints << t_i << endl;
        write2DVectorToFile(fpoints, state_kp);
        return true;
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
};

std::vector<vector<double>> make_web(std::vector<std::vector<double>> matrix, double h_x, double h_y){
    std::vector<std::vector<double>> new_matrix;


    std::vector<double> x_coord;
    x_coord.push_back(0);
    for (int i = 0; i <= (matrix[0]).size()-1; i++){
        x_coord.push_back(h_x * i);
    }
    new_matrix.push_back(x_coord);


    for (int i = 0; i < matrix.size(); i++) {
        std::vector<double> vec;
        vec.push_back(h_y * i);
        for (int j = 0; j < matrix[0].size(); j++) {
            vec.push_back(matrix[i][j]);
        }
        new_matrix.push_back(vec);
    }

    return new_matrix;
}

bool LongTransScheme_for_tables(const PDEProblem &problem, const string &filename, const double & EPS) {

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
        double hx = problem.hx;
        double hy = problem.hy;
        double tau = problem.tau;
        //double tau = std::sqrt(hx*hx + hy*hy)/std::sqrt(1./(problem.X*problem.X) + 1./(problem.Y*problem.Y))/M_PI;
        std::cout<<tau<<std::endl;
        double half_tau = tau/2;

        std::vector<std::vector<double>> state_kp(state_k);
        std::vector<std::vector<double>> state_kph(state_k);

        // Utility vectors
        std::vector<double> Axs(problem.num_x_steps+1, 0.);
        std::vector<double> Bxs(problem.num_x_steps+1, 0.);
        std::vector<double> Cxs(problem.num_x_steps+1, 0.);
        std::vector<double> Dxs(problem.num_x_steps+1, 0.);
        std::vector<double> Ays(problem.num_x_steps+1, 0.);
        std::vector<double> Bys(problem.num_y_steps+1, 0.);
        std::vector<double> Cys(problem.num_y_steps+1, 0.);
        std::vector<double> Dys(problem.num_y_steps+1, 0.);

        // Задаём внешние силы (если нету, то ноль будем прибавлять в коэфах)
        std::function<double(std::vector<double>)> f = ([&] (const std::vector<double>& point) {return  0.;});
        if(problem.extForcesFunction_isSet){
            f = ([&] (const std::vector<double>& point) {return  -problem.extForcesFunction(point);});
        }

        // Максмальное число итераций
        int max_allowed_its = 5000;
        int cur_its = 0;
        do{
            // слой ,,половинный'' первый (k+1/2)
            t_i += half_tau;
            //std::cout << "t_i = " << t_i  << "; cur_its = " << cur_its <<std::endl;
            x_i = problem.x0;
            y_i = problem.y0;

            state_k = state_kp;
            state_kph = state_k;

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
                    Axs[0] = 0.;
                    Bxs[0] = -1./*-(1./tau + 1./(hx*hx))*/;
                    Cxs[0] = -1. / (hx * hx / tau + 1.)/*-1./(hx*hx)*/;
                    Dxs[0] = -(state_k[0][i2] * hx / tau + problem.neymanBoundaryFunc_West({x_i,y_i})) / (hx / tau + 1. / hx)/*-( (1./(2.*hy*hy)) * (state_k[0][i2+1] - 2.*state_k[0][i2] + state_k[0][i2-1]) + (1./tau)*state_k[0][i2] + (1./hx) * problem.neymanBoundaryFunc_West({problem.x0,y_i}) + (0.5)*f({problem.x0,y_i}))*/;
                }

                // Г.У. Восток
                if (problem.dirichletBoundaryFunc_East_isSet) {
                    Axs[problem.num_x_steps] = 0.;
                    Bxs[problem.num_x_steps] = 1.;
                    Cxs[problem.num_x_steps] = 0.;
                    Dxs[problem.num_x_steps] = problem.dirichletBoundaryFunc_East({problem.X, y_i});
                } else if (problem.neymanBoundaryFunc_East_isSet) {
                    /* аппроксимация второго рода*/
                    Axs[problem.num_x_steps] = -1. / (hx * hx / tau + 1.)/*-1./(hx*hx)*/;
                    Bxs[problem.num_x_steps] = -1./*-(1./tau + 1./(hx*hx))*/;
                    Cxs[problem.num_x_steps] = 0.;
                    Dxs[problem.num_x_steps] = -(state_k[problem.num_x_steps][i2] * hx / tau + problem.neymanBoundaryFunc_East({x_i,y_i})) / (hx / tau + 1. / hx) /*-((1./(2*hy*hy))*(state_k[problem.num_x_steps][i2+1] - 2.*state_k[problem.num_x_steps][i2] + state_k[problem.num_x_steps][i2-1]) + (1./hx)*problem.neymanBoundaryFunc_East({problem.X, y_i}) + 1/tau*state_k[problem.num_x_steps][i2] + (1./4.)*(f({problem.X,y_i})+f({problem.X-hx/2.,y_i})))*/;
                }

                for(int i1 = 1; i1 < problem.num_x_steps; ++i1){
                    x_i += hx;
                    Axs[i1] = 1/(hx*hx);
                    Bxs[i1] = 2*( 1/(hx*hx) + 1/tau );
                    Cxs[i1] = 1/(hx*hx);
                    Dxs[i1] = F(state_k[i1][i2], state_k[i1][i2-1],state_k[i1][i2+1], x_i, y_i, problem);
                }

                std::vector<double> x_line_sol = TridiagonalMatrixAlgorithm(Axs, Bxs, Cxs, Dxs);
                for(int i1 = 0; i1<=problem.num_x_steps; ++i1){
                    state_kph[i1][i2] = x_line_sol[i1];
                    //out(state_k[i1]);
                }
                state_kp[0][i2] = state_kph[0][i2];
                state_kp[problem.num_x_steps][i2] = state_kph[problem.num_x_steps][i2];
                //std::cout<<"tridig sol: "<<std::endl;
                //out(x_line_sol);
                //std::cout<<"---"<<std::endl;
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
                    Ays[problem.num_y_steps] = -1. / (hy * hy / tau + 1)/*-1./(hy*hy)*/;
                    Bys[problem.num_y_steps] = -1./*-(1./tau + 1./(hy*hy))*/;
                    Cys[problem.num_y_steps] = 0.;
                    Dys[problem.num_y_steps] = -(state_kph[i1][problem.num_y_steps] * hy / tau + problem.neymanBoundaryFunc_North({x_i, y_i})) / (hy / tau + 1. / hy)/*-((1./(2.*hx*hx))*(state_kph[i1+1][problem.num_y_steps] - 2.*state_kph[i1][problem.num_y_steps] + state_kph[i1-1][problem.num_y_steps]) + (1./hy)*problem.neymanBoundaryFunc_North({x_i, problem.Y}) + 1./tau*state_kph[i1][problem.num_y_steps] + (0.5)*f({x_i,problem.Y}))*/;
                }
                // Г.У. Юг
                if (problem.dirichletBoundaryFunc_South_isSet) {
                    Ays[0] = 0.;
                    Bys[0] = 1.;
                    Cys[0] = 0.;
                    Dys[0] = problem.dirichletBoundaryFunc_South({x_i, problem.y0});
                } else if (problem.neymanBoundaryFunc_South_isSet) {
                    /* аппроксимация второго рода*/
                    Ays[0] = 0.;
                    Bys[0] = -1./*-(1./tau + 1./(hy*hy))*/;
                    Cys[0] = -1. / (hy * hy / tau + 1.)/*-1./(hy*hy)*/;
                    Dys[0] = -(state_kph[i1][0] * hy / tau + problem.neymanBoundaryFunc_South({x_i, y_i})) / (hy / tau + 1. / hy)/*-( (1./(2.*hx*hx)) * (state_kp[i1+1][0] - 2*state_kp[i1][0] + state_kph[i1-1][0]) + (1./tau)*state_kp[i1][0] + 1./tau*state_kph[i1][0] + (1./hy) * problem.neymanBoundaryFunc_South({x_i,problem.y0}) + (0.5)*f({x_i,problem.y0}))*/;
                }

                for(int i2 = 1; i2 < problem.num_y_steps; ++i2){
                    Ays[i2] = 1/(hy*hy);
                    Bys[i2] = 2*( 1/(hy*hy) + 1/tau );
                    Cys[i2] = 1/(hy*hy);
                    Dys[i2] = F_hat(state_kph[i1][i2], state_kph[i1-1][i2],state_kph[i1+1][i2], x_i, y_i, problem);
                }

                std::vector<double> y_line_sol = TridiagonalMatrixAlgorithm(Ays, Bys, Cys, Dys);
                for(int i2 = 0; i2<=problem.num_y_steps; ++i2) {
                    state_kp[i1][i2] = y_line_sol[i2];
                }

                //std::cout<<"tridig sol: "<<std::endl;
                //out(y_line_sol);
                //std::cout<<"---"<<std::endl;
            }
            x_i = problem.x0;
            y_i = problem.y0;
//            fpoints << t_i << endl;
//           write2DVectorToFile(fpoints, state_kp);
            //state_k.swap(state_kp);
            ++cur_its;
        }while(Error_calc(hx*hy*EPS, state_k, state_kph, state_kp) && cur_its < max_allowed_its/*norm1(state_k + (-1)*state_kp) > 1e-17 && cur_its < max_allowed_its*/ );

        state_k = make_web(state_kp, hx, hy);
        write2DVectorToFile(fpoints, state_k);

        return true;
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
};
