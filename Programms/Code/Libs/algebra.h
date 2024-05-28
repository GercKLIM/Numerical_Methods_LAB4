////
//// Объявление функций и переопределений для std::vector
//// для возможности абстракции в математические векторы, матрицы и другой доп. функционал
////
//
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
//
//
///* *** Начальные функции для испорта/экспорта данных *** */
//
//
//
///* Функция импорта матрицы из текстового файла*/
//template <typename T>
//std::vector<std::vector<T>> importSLAU(const std::string& filename);
//
//
///* Функция вывода матрицы на экран */
//template <typename T>
//void print(const std::vector<std::vector<T>>& matrix);
//
//
///* Функция вывода вектора на экран */
//template <typename T>
//void print(const std::vector<T>& vec);
//
//
///* Функция вывода обрезанного вектора на экран */
//template <typename T>
//void print_short(const std::vector<T>& vec, const int& n);
//
//
///* Функция, которая красиво выводит вектор*/
//template<typename T>
//void print_vec(const std::vector<T>& vec);
//
//
///* Функция вывода разделительной линии на экран */
//void printline(const int& n);
//
//
///* Функция для получения матрицы из СЛАУ */
//template <typename T>
//std::vector<std::vector<T>> SLAU_to_matrix(const std::vector<std::vector<T>>& SLAU);
//
//
///* Функция для получения векторая из СЛАУ */
//template <typename T>
//std::vector<T> SLAU_to_vec(const std::vector<std::vector<T>>& SLAU);
//
//
//
///* *** Функции математики векторов *** */
//
//
//
///* Операция cложения векторов */
//template <typename T>
//std::vector<T> operator+(const std::vector<T>& vec1, const  std::vector<T>& vec2);
//
//
///* Операция вычитания векторов */
//template <typename T>
//std::vector<T> operator-(const std::vector<T>& vec1, const std::vector<T>& vec2);
//
//
///* Операция почленного умножения векторов */
//template <typename T>
//std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2);
//
//
///* Операция умножения вектора на число */
//template <typename T>
//std::vector<T> operator*(const T& c, const std::vector<T>& vec2);
//
//
//template <typename T>
//std::vector<T> operator*(const std::vector<T>& vec2, const T& c);
//
///* Операция деления вектора на число */
//template<typename T>
//std::vector<T> operator/(const std::vector<T>& vec, const T& c);
//
//
///* Операция почленного деления векторов */
//template <typename T>
//std::vector<T> operator/(const std::vector<T>& vec1, const std::vector<T>& vec2);
//
//
//// Определение оператора отрицания для матрицы
//template <typename T>
//std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& matrix);
//
//
///* Функция для скалярного умножения векторов */
//template <typename T>
//T dot(const std::vector<T>& vec1, const std::vector<T>& vec2);
//
//
///* Функция для нормы вектора */
//template <typename T>
//T norm(const std::vector<T>& vec, const int& p = 2);
//
//
///* Функция, которая возращает матрицу комбинаций элементов вектора */
//template<typename T>
//std::vector<std::vector<T>> generateCombinations(const std::vector<T>& vec);
//
//
///* Функция, возвращает вектор модулей */
//template<typename T>
//std::vector<T> vec_abs(const std::vector<T>& vec);
//
//
///* Функция, возращающая сумму элементов вектора */
//template<typename T>
//T sum(const std::vector<T>& vec);
//
//
//
//
///* *** Функции математики матриц *** */
//
//
//
//
///* Матричное умножение */
//template <typename T>
//std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);
//
//
///* Функция поэлементного сложения матриц */
//template <typename T>
//std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);
//
//
///* Функция поэлементного вычитания матриц */
//template <typename T>
//std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);
//
//
///* Функция для умножения матрицы на вектор */
//template <typename T>
//std::vector<T> operator*(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec);
//
//
///* Функция для транспонирования матрицы */
//template <typename T>
//std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& A);
//
//
///* Функция для создания единичной матрицы размера n x n */
//template <typename T>
//std::vector<std::vector<T>> create_identity_matrix(const int& n);
//
//
///* Функция для поэлементного умножения матриц */
//template <typename T>
//std::vector<std::vector<T>> Multyply(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);
//
//
///* Функция округления чисел в матрицах */
//template <typename T>
//std::vector<std::vector<T>> Matrix_round(const std::vector<std::vector<T>>& A, const double& eps);
//
//
///* Функция для вычисления нормы матрицы */
//template <typename T>
//T norm(const std::vector<std::vector<T>>& matrix, const int& p = 2);
//
//
///* Функция для вычисления числа обусловленности матрицы c нормой 1*/
//template <typename T>
//T cond_1(const std::vector<std::vector<T>>& matrix);
//
//
///* Функция для вычисления числа обусловленности матрицы c нормой 2*/
//template <typename T>
//T cond_2(const std::vector<std::vector<T>>& matrix);
//
//
///* Функция для вычисления числа обусловленности матрицы c нормой oo*/
//template <typename T>
//T cond_oo(const std::vector<std::vector<T>>& matrix);
//
//
///* Функция поворота матрицы вправо */
//template <typename T>
//std::vector<std::vector<T>> RotateRight(const std::vector<std::vector<T>>& A);
//
//
///* Функция поворота матрицы влево */
//template <typename T>
//std::vector<std::vector<T>> RotateLeft(const std::vector<std::vector<T>>& A);
//
//
//// Функция для обратной матрицы с проверкой на вырожденность c определенной точностью
//template <typename T>
//std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A, const T& eps);
//
//
//// Функция для обратной матрицы с проверкой на вырожденность c определенной точностью
//template <typename T>
//std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A);
//
//
//// Функция обрезки матрицы снизу и справа
//template <typename T>
//std::vector<std::vector<T>> crop_matrix(const std::vector<std::vector<T>>& A, const int& k);
//
//
///* Функция, вычисляющая определитель матрицы 4х4 */
//template <typename T>
//double det(const std::vector<std::vector<T>>& matrix);
//
//
///* Функция, сортирующая вектор */
//template< typename T>
//std::vector<T> sorted(const std::vector<T>& vec_not_sort);
//
//
///* Функция, возращающая максимальный элемент вектора */
//template<typename T>
//T vec_max(const std::vector<T>& vec);
//
//
///* Функция, вычисляющая норму разности векторов */
//double sqr(std::vector<double> vec1, std::vector<double> vec2);
//
//
///* Функция, численно вычисляющая произвоную в точке point по i переменной */
//double Differential(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const int& i, const double& eps);
//
//
///* Функция, вычисляющая градиент функции в точке point */
//std::vector<double> Gradient(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const double& eps);
//
//
///* Функция для сдвига вектора на n элементов */
//template<typename T>
//std::vector<T> shift(const std::vector<T>& vec, int n);

/*<<ќперации линейной агебры (матрицами и веторы)>>
 @author: Ўаманов »ван
 @author:  лимов ќлег
 2021-2024*/
#pragma once
#include<vector>
#include<iomanip>
#include<iostream>
using namespace std;

// сложение матриц
template<typename LT>
vector<vector<LT>> mat_sum(const vector<vector<LT>>& matr1, const vector<vector<LT>>& matr2)
{
    vector<vector <LT>> product(matr1);
    for (int i = 0; i < matr1.size(); ++i)
        for (int j = 0; j < matr2.size(); ++j)
            product[i][j] += matr2[i][j];
    return product;
}
template<typename LT>
vector<vector<LT>> operator + (const vector<vector<LT>> matr1, const vector<vector<LT>> matr2)
{
    return mat_sum(matr1, matr2);
}
//умножение матрицы на матрицу
template<typename LT>
vector<vector<LT>> mat_mult(const vector<vector<LT>> matr1, const vector<vector<LT>> matr2)
{
    vector<vector<LT>> product;
    int n = matr1[0].size();
    product.push_back({});
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            LT temp_sum = 0;
            for (int k = 0; k < n; ++k)
                temp_sum += matr1[i][k] * matr2[k][j];
            product[i].push_back(temp_sum);
        }
        product.push_back({});
    }
    return product;
}

template<typename LT>
vector<vector<LT>> operator * (const vector<vector<LT>> matr1, const vector<vector<LT>> matr2)
{
    return mat_mult(matr1, matr2);
}

//умножение матрицы на число
template<typename LT, typename LT2>
vector<vector<LT>> operator * (const vector<vector<LT>>& matr, const LT2 lambda)
{
    vector<vector<LT>> result(matr);
    for (int i = 0; i < matr.size(); ++i)
        for (int j = 0; j < matr[0].size(); ++j)
            result[i][j] *= lambda;
    return result;
}
template<typename LT, typename LT2>
vector<vector<LT>> operator * (const LT2 lambda, const vector<vector<LT>>& matr)
{
    return matr * lambda;
}
//умножение матрицы на вектор
template<typename LT>
vector<LT> operator * (const vector<vector<LT>>& matr, const vector<LT> vec)
{
    vector<LT> result;
    for (int i = 0; i < vec.size(); ++i)
    {
        LT tmp_val = 0;
        for (int j = 0; j < vec.size(); ++j)
        {
            tmp_val += matr[i][j] * vec[j];
        }
        result.push_back(tmp_val);
    }
    return result;
}
//умножение вектора на число
template<typename LT, typename LT2>
vector<LT> operator * (const vector<LT>& vec, const LT2& lambda)
{
    vector<LT> result;
    for (auto el : vec)
    {
        result.push_back(el * lambda);
    }
    return result;
}
template<typename LT, typename LT2>
vector<LT> operator * (const LT2& lambda, const vector<LT>& vec)
{
    return vec * lambda;
}

//скал¤рное произведение
template <typename LT>
LT scmult(vector<LT> v1, vector<LT> v2)
{
    LT result = 0;
    for (auto v = v1.begin(), p = v2.begin(); v != v1.end() && p != v2.end(); v++, p++)
    {
        result += *v * *p;
    }
    return result;
}
template<typename LT>
LT operator * (const vector<LT>& v1, const vector<LT>& v2)
{
    return scmult(v1, v2);
}
//норма
template <typename LT>
LT vec_norm(vector<LT> v1)
{
    return sqrt(scmult<LT>(v1, v1));
}
//cумма
template <typename LT>
vector<LT> vsum(vector<LT>& v1, vector<LT>& v2)
{
    vector<LT> result;
    auto v = v1.begin();
    auto p = v2.begin();
    for (; v != v1.end() && p != v2.end(); v++, p++)
    {
        result.push_back(*v + *p);
    }
    for (; v != v1.end(); v++)
    {
        result.push_back(*v);
    }
    for (; p != v2.end(); p++)
    {
        result.push_back(*p);
    }
    return result;
}
//разность
template <typename LT>
vector<LT> vdiff(vector<LT>& v1, vector<LT>& v2)
{
    vector<LT> result;
    auto v = v1.begin();
    auto p = v2.begin();
    for (; v != v1.end() && p != v2.end(); v++, p++)
    {
        result.push_back(*v - *p);
    }
    for (; v != v1.end(); v++)
    {
        result.push_back(*v);
    }
    for (; p != v2.end(); p++)
    {
        result.push_back(-*p);
    }
    return result;
}

//сложение векторов
template<typename LT>
vector<LT> operator + (vector<LT> r_vec, vector<LT> l_vec)
{
    return vsum(r_vec, l_vec);
}
//вычитание векторов
template<typename LT>
vector<LT> operator - (vector<LT> r_vec, vector<LT> l_vec)
{
    return vdiff(r_vec, l_vec);
}

//прибавление единичной*lambda  к матрице
template<typename LT, typename LT2>
vector<vector<LT>> operator + (vector<vector<LT>> matr, LT2 lambda)
{
    vector<vector<LT>> result(matr);
    for (int i = 0; i < matr[0].size(); ++i)
        result[i][i] += lambda;
    return result;
}
template<typename LT, typename LT2>
vector<vector<LT>> operator + (LT2 lambda, vector<vector<LT>> matr)
{
    return matr + lambda;
}
//Ќќ–ћџ ћј“–»÷
template<typename LT>
LT norm1(vector<vector<LT>> A)
{
    int n = A[0].size();
    LT max_val = 0;
    for (int j = 0; j < n; ++j)
    {
        LT temp_sum = 0;
        for (int i = 0; i < n; ++i)
        {
            temp_sum += fabs(A[i][j]);
        }
        if (temp_sum > max_val)
            max_val = temp_sum;
    }
    return max_val;
}

template<typename LT>
LT norm_inf(vector<vector<LT>> A)
{
    int n = A[0].size();
    LT max_val = 0;
    for (int i = 0; i < n; ++i)
    {
        LT temp_sum = 0;
        for (int j = 0; j < n; ++j)
        {
            temp_sum += fabs(A[i][j]);
        }
        if (temp_sum > max_val)
            max_val = temp_sum;
    }
    return max_val;
}
template<typename LT>
LT vec_norm_inf(vector<LT> vec)
{
    LT max_el = 0;
    for (auto el : vec)
        if (max_el < fabs(el))
            max_el = fabs(el);
    return max_el;
}

template<typename LT>
LT vec_norm1(vector<LT> vec)
{
    LT res = 0;
    for (auto el : vec)
        res += fabs(el);
    return res;
}

//печать вектора
template<typename LT>
void out(vector<LT> vec)
{
    int n = vec.size();
    for (int i = 0; i < n; ++i)
    {
        cout << fixed << setprecision(6) << setw(8) << setfill(' ') << vec[i] << "  ";
    }
    cout << endl;

}

template<typename LT>
LT det(vector<vector<LT>> matrix)
{
    int dim_size = matrix.size();
    if (dim_size == matrix[0].size())
    {
        if (dim_size == 2)
        {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        else if (dim_size == 3)
        {
            return	matrix[0][0] * matrix[1][1] * matrix[2][2] +
                      matrix[2][0] * matrix[0][1] * matrix[1][2] +
                      matrix[0][2] * matrix[1][0] * matrix[2][1] -
                      matrix[2][0] * matrix[1][1] * matrix[0][2] -
                      matrix[1][0] * matrix[0][1] * matrix[2][2] -
                      matrix[0][0] * matrix[2][1] * matrix[1][2];
        }
    }
    return;
}

//транспонирование матрицы
template<typename LT>
vector<vector<LT>> transpose(vector<vector<LT>> matr)
{
    int n = matr[0].size();
    vector<vector<LT>> product(matr);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (i != j)
            {
                product[i][j] = matr[j][i];
                product[j][i] = matr[i][j];
            }
    return product;
}