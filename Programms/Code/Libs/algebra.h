//
// Объявление функций и переопределений для std::vector
// для возможности абстракции в математические векторы, матрицы и другой доп. функционал
//

#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>


/* *** Начальные функции для испорта/экспорта данных *** */



/* Функция импорта матрицы из текстового файла*/
template <typename T>
std::vector<std::vector<T>> importSLAU(const std::string& filename);


/* Функция вывода матрицы на экран */
template <typename T>
void print(const std::vector<std::vector<T>>& matrix);


/* Функция вывода вектора на экран */
template <typename T>
void print(const std::vector<T>& vec);


/* Функция вывода обрезанного вектора на экран */
template <typename T>
void print_short(const std::vector<T>& vec, const int& n);


/* Функция, которая красиво выводит вектор*/
template<typename T>
void print_vec(const std::vector<T>& vec);


/* Функция вывода разделительной линии на экран */
void printline(const int& n);


/* Функция для получения матрицы из СЛАУ */
template <typename T>
std::vector<std::vector<T>> SLAU_to_matrix(const std::vector<std::vector<T>>& SLAU);


/* Функция для получения векторая из СЛАУ */
template <typename T>
std::vector<T> SLAU_to_vec(const std::vector<std::vector<T>>& SLAU);



/* *** Функции математики векторов *** */



/* Операция cложения векторов */
template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const  std::vector<T>& vec2);


/* Операция вычитания векторов */
template <typename T>
std::vector<T> operator-(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Операция почленного умножения векторов */
template <typename T>
std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Операция умножения вектора на число */
template <typename T>
std::vector<T> operator*(const T& c, const std::vector<T>& vec2);


template <typename T>
std::vector<T> operator*(const std::vector<T>& vec2, const T& c);

/* Операция деления вектора на число */
template<typename T>
std::vector<T> operator/(const std::vector<T>& vec, const T& c);


/* Операция почленного деления векторов */
template <typename T>
std::vector<T> operator/(const std::vector<T>& vec1, const std::vector<T>& vec2);


// Определение оператора отрицания для матрицы
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& matrix);


/* Функция для скалярного умножения векторов */
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Функция для нормы вектора */
template <typename T>
T norm(const std::vector<T>& vec, const int& p = 2);


/* Функция, которая возращает матрицу комбинаций элементов вектора */
template<typename T>
std::vector<std::vector<T>> generateCombinations(const std::vector<T>& vec);


/* Функция, возвращает вектор модулей */
template<typename T>
std::vector<T> vec_abs(const std::vector<T>& vec);


/* Функция, возращающая сумму элементов вектора */
template<typename T>
T sum(const std::vector<T>& vec);




/* *** Функции математики матриц *** */




/* Матричное умножение */
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция поэлементного сложения матриц */
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция поэлементного вычитания матриц */
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция для умножения матрицы на вектор */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec);


/* Функция для транспонирования матрицы */
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& A);


/* Функция для создания единичной матрицы размера n x n */
template <typename T>
std::vector<std::vector<T>> create_identity_matrix(const int& n);


/* Функция для поэлементного умножения матриц */
template <typename T>
std::vector<std::vector<T>> Multyply(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция округления чисел в матрицах */
template <typename T>
std::vector<std::vector<T>> Matrix_round(const std::vector<std::vector<T>>& A, const double& eps);


/* Функция для вычисления нормы матрицы */
template <typename T>
T norm(const std::vector<std::vector<T>>& matrix, const int& p = 2);


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(const std::vector<std::vector<T>>& matrix);


/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(const std::vector<std::vector<T>>& matrix);


/* Функция для вычисления числа обусловленности матрицы c нормой oo*/
template <typename T>
T cond_oo(const std::vector<std::vector<T>>& matrix);


/* Функция поворота матрицы вправо */
template <typename T>
std::vector<std::vector<T>> RotateRight(const std::vector<std::vector<T>>& A);


/* Функция поворота матрицы влево */
template <typename T>
std::vector<std::vector<T>> RotateLeft(const std::vector<std::vector<T>>& A);


// Функция для обратной матрицы с проверкой на вырожденность c определенной точностью
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A, const T& eps);


// Функция для обратной матрицы с проверкой на вырожденность c определенной точностью
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A);


// Функция обрезки матрицы снизу и справа
template <typename T>
std::vector<std::vector<T>> crop_matrix(const std::vector<std::vector<T>>& A, const int& k);


/* Функция, вычисляющая определитель матрицы 4х4 */
template <typename T>
double det(const std::vector<std::vector<T>>& matrix);


/* Функция, сортирующая вектор */
template< typename T>
std::vector<T> sorted(const std::vector<T>& vec_not_sort);


/* Функция, возращающая максимальный элемент вектора */
template<typename T>
T vec_max(const std::vector<T>& vec);


/* Функция, вычисляющая норму разности векторов */
double sqr(std::vector<double> vec1, std::vector<double> vec2);


/* Функция, численно вычисляющая произвоную в точке point по i переменной */
double Differential(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const int& i, const double& eps);


/* Функция, вычисляющая градиент функции в точке point */
std::vector<double> Gradient(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const double& eps);


/* Функция для сдвига вектора на n элементов */
template<typename T>
std::vector<T> shift(const std::vector<T>& vec, int n);