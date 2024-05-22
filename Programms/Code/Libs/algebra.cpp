//
// Реализация функций и переопределений для std::vector
// для возможности абстракции в математические векторы, матрицы и т.д.
//

#include "algebra.h"



/* *** Начальные функции для испорта/экспорта данных *** */

/* Функция импорта чисел из файла */
std::vector<double> ImportData(const std::string& filename) {
    std::vector<double> numbers;

    // Открытие файла для чтения
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return numbers;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Преобразование строки в число типа double и добавление его в вектор
        double num = std::atof(line.c_str());
        numbers.push_back(num);
    }

    // Закрытие файла
    file.close();

    return numbers;
}

/* Функция импорта матрицы из текстового файла*/
template <typename T>
std::vector<std::vector<T>> importSLAU(const std::string& filename) {
    std::vector<std::vector<T>> matrix;
    std::vector<T> vec;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << "Error: not open file \n" << std::endl;
        exit(1);
    }

    int size;
    file >> size;

    matrix.resize(size, std::vector<T>(size+1));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size+1; ++j) {
            T value;
            if (file >> value) {
                matrix[i][j] = value;
            }
        }
    }

    file.close();
    return matrix;
};


/* Функция вывода матрицы на экран */
template <typename T>
void print(const std::vector<std::vector<T>>& matrix) {
    for (std::vector<T> row : matrix) {
        for (T value: row) {
            std::cout << value << ' ';
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}


/* Функция вывода вектора на экран */
template <typename T>
void print(const std::vector<T>& vec) {
    for (T value : vec) {
        std::cout << value << ' ';
    }
    std::cout << std::endl;
}


/* Функция вывода обрезанного вектора на экран */
template <typename T>
void print_short(const std::vector<T>& vec, const int& n){

    for (int i = 0; i < n; ++i){
        std::cout << vec[i] << ' ';
    }
    std::cout << "..." << std::endl;
}


/* Функция, которая красиво выводит вектор*/
template<typename T>
void print_vec(const std::vector<T>& vec){
    std::cout << "(" << vec[0];
    for (int i = 1; i < vec.size(); i++){
        std::cout << ", " << vec[i];
    }
    std::cout << ")" << std::endl;
}


/* Функция вывода разделительной линии на экран */
void printline(const int& n){
    for (int i = 0; i < n; i ++){
        std::cout << "-";
    }
    std::cout << std::endl;
}


/* Функция для получения матрицы из СЛАУ */
template <typename T>
std::vector<std::vector<T>> SLAU_to_matrix(const std::vector<std::vector<T>>& SLAU){
    std::vector<std::vector<T>> matrix;
    matrix.resize(SLAU.size(), std::vector<T>(SLAU.size()));

    for (int i = 0; i < SLAU.size(); i++) {
        for (int j = 0; j < SLAU.size(); j++) {
            matrix[i][j] = SLAU[i][j];
        }
    }
    return matrix;
}


/* Функция для получения вектора из СЛАУ */
template <typename T>
std::vector<T> SLAU_to_vec(const std::vector<std::vector<T>>& SLAU){
    int s = SLAU.size();
    std::vector<T> vec(s);

    for (int i = 0; i < SLAU.size(); i++) {
        vec[i] = SLAU[i][s];
    }
    return vec;
}



/* *** Функции математики векторов *** */



/* Функция для сложения векторов */
template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2){
    std::vector<T> pert_vec = vec1;
    for (int i = 0; i < vec1.size(); i++) {
        pert_vec[i] += vec2[i];
    }
    return pert_vec;
}


/* Функция вычитания векторов */
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
    // Проверка на возможность умножения
    if (a.size() != b.size()) {
        std::cout << "Error: size a != size b in substraction vectors." << std::endl;
        exit(1);
    }
    // Создание результирующего вектора
    std:: vector<T> result(a.size(), 0);

    // Умножение матрицы на вектор
    for (int i = 0; i < a.size(); ++i) {
        result[i] += a[i] - b[i];
    }
    return result;

}


/* Операция почленного умножения векторов */
template <typename T>
std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in operator*." << std::endl;
        exit(1);
    }
    std::vector<T> result(vec1.size(), 0);
    for (int i = 0; i < vec1.size(); i++){
        result[i] = vec1[i] * vec2[i];
    }
    return result;
}


/* Операция умножения вектора на число */
template <typename T>
std::vector<T> operator*(const T& c, const std::vector<T>& vec){
    std::vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, const T& c){
    std::vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}

/* Операция деления вектора на число */
template<typename T>
std::vector<T> operator/(const std::vector<T>& vec, const T& c) {
    std::vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++) {
        result[i] = vec[i] / c;
    }
    return result;
}


/* Операция почленного деления векторов */
template <typename T>
std::vector<T> operator/(const std::vector<T>& vec1, const std::vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in operator*." << std::endl;
        exit(1);
    }
    std::vector<T> result(vec1.size(), 0);
    for (int i = 0; i < vec1.size(); i++){
        result[i] = vec1[i] / vec2[i];
    }
    return result;
}


/* Функция для скалярного умножения векторов */
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in operator*." << std::endl;
        exit(1);
    }
    T result;
    for (int i = 0; i < vec1.size(); i++){
        result += vec1[i] * vec2[i];
    }
    return result;
}


/* Функция для нормы вектора */
template <typename T>
T norm(const std::vector<T>& vec, const int& p){
    if (vec.empty()) {
        std::cerr << "Error: Empty vector in norm() \n";
        exit(1);
    }

    T result = 0.0;

    // Вычисление нормы
    if (p == 0) {
        // Норма oo
        for (const auto& element : vec) {
            T absElement = abs(element);
            if (absElement > result) {
                result = absElement;
            }
        }
    } else {
        // Общий случай для норм L1, L2 и т.д.
        for (const auto& element : vec) {
            result += pow(abs(element), p);
        }

        result = pow(result, 1.0 / p);
    }

    return result;
}


/* Функция, которая возращает матрицу комбинаций элементов вектора */
template<typename T>
std::vector<std::vector<T>> generateCombinations(const std::vector<T>& vec) {
    int n = vec.size();

    // Вектор для хранения всех комбинаций
    std::vector<std::vector<T>> combinations;

    // Внешний цикл по всем возможным комбинациям
    for (int i = 0; i < (1 << n); ++i) {
        std::vector<T> current(n);

        // Внутренний цикл для каждой позиции вектора
        for (int j = 0; j < n; ++j) {
            current[j] = (i & (1 << j)) ? vec[j] : -vec[j];
        }

        // Добавить текущую комбинацию в вектор
        combinations.push_back(current);
    }

    return combinations;
}


/* Функция, возвращает вектор модулей */
template<typename T>
std::vector<T> vec_abs(const std::vector<T>& vec){
    for (int i = 0; i < vec.size(); i++){
        vec[i] = fabs(vec[i]);
    }
    return vec;
}


/* Функция, возращающая сумму элементов вектора */
template<typename T>
T sum(const std::vector<T>& vec){
    T sum = 0;
    for (int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum;
}



/* *** Функции математики матриц *** */



/* Операция для умножения матрицы на число */
//template <typename T>
//vector<vector<T>> operator*(const vector<vector<T>>& A, const T& scalar){
//    // Создание результирующей матрицы с теми же размерами
//    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));
//
//    // Умножение каждого элемента матрицы на число
//    for (size_t i = 0; i < A.size(); ++i) {
//        for (size_t j = 0; j < A[0].size(); ++j) {
//            result[i][j] = A[i][j] * scalar;
//        }
//    }
//
//    return result;
//}


/* Операция для умножения  числа на матрицу */
//template <typename T>
//vector<vector<T>> operator*(const T& scalar, const vector<vector<T>>& A){
//    // Создание результирующей матрицы с теми же размерами
//    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));
//
//    // Умножение каждого элемента матрицы на число
//    for (size_t i = 0; i < A.size(); ++i) {
//        for (size_t j = 0; j < A[0].size(); ++j) {
//            result[i][j] = A[i][j] * scalar;
//        }
//    }
//
//    return result;
//}


/* Операция поэлементного сложения матриц */
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    // Проверка на совпадение размеров матриц
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        std::cout << "Error: size A != size B in addition matrix." << std::endl;
        exit(1);
    }

    // Создание результирующей матрицы с теми же размерами
    std::vector<std::vector<T>> result(A.size(), std::vector<T>(A[0].size(), 0));

    // Поэлементное сложение
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}


/* Операция поэлементного вычитания матриц */
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    // Проверка на совпадение размеров матриц
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        std::cout << "Error: size A != size B in substraction matrix." << std::endl;
        exit(1);
    }

    // Создание результирующей матрицы с теми же размерами
    std::vector<std::vector<T>> result(A.size(), std::vector<T>(A[0].size(), 0));

    // Поэлементное сложение
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}


/* Операция умножения матрицы на вектор */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec) {
    // Проверка на возможность умножения
    if (matrix[0].size() != vec.size()) {
        std::cout << "Error: size A != size b in multiply Matrix By Vector." << std::endl;
        exit(1);
    }
    // Создание результирующего вектора
    std::vector<T> result(matrix.size(), 0);

    // Умножение матрицы на вектор
    for (int i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}


/* Матричное умножение */
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    int m = A.size();    // Количество строк в матрице A
    int n = A[0].size(); // Количество столбцов в матрице A
    int p = B[0].size(); // Количество столбцов в матрице B

    if (n != B.size()) {
        std::cout << "Error: impossible multiply matrix" << std::endl;
        exit(1);
    }

    std::vector<std::vector<T>> result(m, std::vector<T>(p, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}


// Определение оператора отрицания для матрицы
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = -matrix[i][j];
        }
    }
    return result;
}


/* Функция для поэлементного умножения матриц */
template <typename T>
std::vector<std::vector<T>> Multyply(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    int m = A.size();    // Количество строк в матрице A
    int n = A[0].size(); // Количество столбцов в матрице A
    int p = B[0].size(); // Количество столбцов в матрице B

    if (n != B.size()) {
        printf("Error: impossible multiply matrix");
        exit(1);
    }

    std::vector<std::vector<T>> result(m, std::vector<T>(p, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            result[i][j] = A[i][j] * B[i][j];
        }
    }
    return result;
}


/* Функция округления чисел в матрицах */
template <typename T>
std::vector<std::vector<T>> Matrix_round(const std::vector<std::vector<T>>& A, const double& eps){
    std::vector<std::vector<T>> roundA = A;
    int size = A.size();

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            roundA[i][j] = (round(A[i][j]) >= 0)? round(abs(A[i][j]) * (1 / eps)) / (1 / eps): -1 * round(abs(A[i][j]) * (1 / eps)) / (1 / eps);
        }
    }
    return roundA;
}


/* Функция для вычисления нормы матрицы */
template <typename T>
T norm(const std::vector<std::vector<T>>& matrix, const int& p) {
    // Проверка на пустую матрицу
    if (matrix.empty() || matrix[0].empty()) {
        std::cout << "Error: Empty matrix in norm()\n" << std::endl;
        exit(1);
    }

    int rows = matrix.size();
    int cols = matrix[0].size();

    T result = 0.0;

    // Вычисление нормы матрицы
    if (p == 0) {
        // Норма матрицы Чебышева (максимальное значение по модулю в строке)
        for (int i = 0; i < rows; ++i) {
            T rowSum = 0.0;
            for (int j = 0; j < cols; ++j) {
                rowSum += abs(matrix[i][j]);
            }
            if (rowSum > result) {
                result = rowSum;
            }
        }
    } else {
        // Общий случай для норм матрицы (Фробениуса и др.)
        for (int j = 0; j < cols; ++j) {
            T colSum = 0.0;
            for (T i = 0; i < rows; ++i) {
                colSum += pow(abs(matrix[i][j]), p);
            }
            result += pow(colSum, 1.0 / p);
        }

        result = pow(result, 1.0 / p);
    }
    return result;
}


/* Функция поворота матрицы вправо */
template <typename T>
std::vector<std::vector<T>> RotateRight(const std::vector<std::vector<T>>& A){

    std::vector<std::vector<T>> A_rotate(A.size(), std::vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[A.size() - 1 - j][i] = A[i][j];
        }
    }

    return A_rotate;

}


/* Функция поворота матрицы влево */
template <typename T>
std::vector<std::vector<T>> RotateLeft(const std::vector<std::vector<T>>& A){

    std::vector<std::vector<T>> A_rotate(A.size(), std::vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[j][A.size() - 1 - i] = A[i][j];
        }
    }

    return A_rotate;
}


// Функция для создания единичной матрицы размера n x n
template <typename T>
std::vector<std::vector<T>> create_identity_matrix(const int& n) {
    std::vector<std::vector<T>> identity(n, std::vector<T>(n, 0));
    for (int i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}


// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A, const T& eps) {
    std::vector<std::vector<T>> E = create_identity_matrix<T>(A.size());
    std::vector<std::vector<T>> E_rotate = RotateLeft(E);
    std::vector<T> e(A.size());
    std::vector<std::vector<T>> X(A.size(), std::vector<T>(A.size(), 0));


    for (int i = 0; i < A.size(); i++){
        e = E_rotate[i];
        X[i] = method_Gaussa(A, e, eps);

    }
    std::vector<std::vector<T>> A_inv = RotateLeft(X);
    return A_inv;
}


// Функция для обратной матрицы с проверкой на вырожденность с максимальной точностью
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A){
    T eps = std::numeric_limits<T>::epsilon();
    return inverseMatrix(A, eps);
}


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(const std::vector<std::vector<T>>& matrix){
    T n_1 = norm_1(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_1(A) = oo");
        return std::numeric_limits<T>::infinity();
    }
    std::vector<std::vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_1(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}


/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(const std::vector<std::vector<T>>& matrix){
    T n_1 = norm_2(matrix);
    if (n_1 == 0) {
        std::cout << "Error: Det(A) = 0  =>  cond_2(A) = oo" << std::endl;
        return std::numeric_limits<T>::infinity();
    }
    std::vector<std::vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_2(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}


/* Функция для вычисления числа обусловленности матрицы с нормой oo*/
template <typename T>
T cond_oo(const std::vector<std::vector<T>>& matrix){
    T n_1 = norm_oo(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_oo(A) = oo");
        return std::numeric_limits<T>::infinity();
    }
    std::vector<std::vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_oo(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}

/* Функция транспонирования матрицы */
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& A) {
    int rows = A.size();
    int cols = A[0].size();
    std::vector<std::vector<T>> result(cols, std::vector<T>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = A[i][j];
        }
    }

    return result;
}


// Функция обрезки матрицы снизу и справа
template <typename T>
std::vector<std::vector<T>> crop_matrix(const std::vector<std::vector<T>>& A, const int& k){

    int n = A.size();
    std::vector<std::vector<T>> A_crop(n - k, std::vector<T>(n - k, 0));
    for (int i = 0; i < (n - k); i++){
        for (int j = 0; j < (n - k); j++){
            A_crop[i][j] = A[i][j];
        }
    }

    return A_crop;
}

/* Функция, вычисляющая определитель матрицы 4х4 */
template <typename T>
double det(const std::vector<std::vector<T>>& matrix) {
    return
            matrix[0][0] * (
                    matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
                    matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +
                    matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1])
            ) -
            matrix[0][1] * (
                    matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
                    matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
                    matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0])
            ) +
            matrix[0][2] * (
                    matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) -
                    matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
                    matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])
            ) -
            matrix[0][3] * (
                    matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) -
                    matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +
                    matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]));
}


/* Функция, сортирующая вектор */
template< typename T>
std::vector<T> sorted(const std::vector<T>& vec_not_sort) {
    std::vector<T> vec(vec_not_sort);
    int n = vec.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (vec[j] > vec[j + 1]) {
                // Обмен элементов, если они не упорядочены
                T temp = vec[j];
                vec[j] = vec[j + 1];
                vec[j + 1] = temp;
            }
        }
    }
    return vec;
}


/* Функция, возращающая максимальный по модулю элемент вектора */
template<typename T>
T vec_max(const std::vector<T>& vec){
    int n = vec.size();
    T max = 0;
    for (int i = 0; i < n; i++) {
        if (abs(vec[i]) > max)
            max = abs(vec[i]);
    }
    return max;
}


/* ### Переопределение потока вывода для vector ### */

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[";
    for (int i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}


/* Функция, вычисляющая норму разности векторов */
double sqr(std::vector<double> vec1, std::vector<double> vec2) {
    int m = vec1.size();
    double sum;
    for (int i = 0; i < m; i++) {
        sum = (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sum;
}


/* Функция, численно вычисляющая произвоную в точке point по i переменной */
double Differential(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const int& i, const double& eps){

    std::vector<double> left_point(point);
    left_point[i] -= eps;
    std::vector<double> right_point(point);
    right_point[i] += eps;

    return (F(right_point)[i] - F(left_point)[i]) / (2 * eps);
}

/* Функция, вычисляющая градиент функции в точке point */
std::vector<double> Gradient(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const double& eps){

    int N = point.size();
    std::vector<double> grad(N, 0);

    for (int i = 0; i < N; i++){
        grad[i] = Differential(F, point, i, eps);
    }
    return grad;
}


/* Функция для сдвига вектора на n элементов */
template<typename T>
std::vector<T> shift(const std::vector<T>& vec, int n) {
    std::vector<T> shiftedVec(vec.size()); // Создаем вектор той же длины
    int size = vec.size();

    // Если сдвиг больше длины вектора, находим остаток от деления
    n = n % size;

    // Перемещаем элементы вправо
    if (n >= 0) {
        for (int i = 0; i < size; ++i) {
            shiftedVec[(i + n) % size] = vec[i];
        }
    }
        // Перемещаем элементы влево
    else {
        for (int i = 0; i < size; ++i) {
            shiftedVec[(size + i + n) % size] = vec[i];
        }
    }

    return shiftedVec;
}

