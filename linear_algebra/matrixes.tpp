#include "matrixes.h"
#include "gauss_method.h"
#include <cmath>
#include <random>

// Конструктор пустой матрицы
template <class T>
matrixes<T>::matrixes(int rows, int cols) : matrixes(nullptr, rows, cols) {}

// Основной конструктор с инициализацией
template <class T>
matrixes<T>::matrixes(T *data, int rows, int cols) : rows(rows), cols(cols) {
    matrix = new T[(rows * cols)];

    if (data == nullptr) {
        for (int i = 0; i < (rows * cols); i++) {
            matrix[i] = T(0.0);
        }
    } else {
        for (int i = 0; i < (rows * cols); i++) {
            matrix[i] = data[i];
        }
    }
}

// Деструктор
template <class T>
matrixes<T>::~matrixes() {
    if (matrix) {
        delete[] matrix;
        matrix = nullptr;
    }
}

// Геттеры для доступа к элементам и размерам
template <class T>
T matrixes<T>::GetIJ(int i, int j) const { return matrix[cols * i + j]; }

template <class T>
T matrixes<T>::Get(int index) const { return matrix[index]; }

template <class T>
T *matrixes<T>::GetMatrix() const { return matrix; }

template <class T>
int matrixes<T>::GetRows() const { return rows; }

template <class T>
int matrixes<T>::GetCols() const { return cols; }

// Заполняет матрицу значениями матрицы Гильберта
template <class T>
void matrixes<T>::MakeGilbert() {
    if (rows != cols) {
        return;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[cols * i + j] = T(1.0) / T(i + j + 1);
        }
    }
}

// Возвращает обратную матрицу
template <class T>
T *matrixes<T>::GetInverseMatrix() const {
    gauss_method<T> mat(matrix, rows);
    mat.TakeReverse();

    T *result = new T[rows * rows];

    for (int i = 0; i < (rows * rows); i++) {
        result[i] = mat.GetMatrix()[i];
    }

    return result;
}

// Создает единичную матрицу
template <class T>
void matrixes<T>::MakeOnes() {
    if (rows != cols) {
        return;
    }

    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = T(0.0);
    }

    for (int i = 0; i < rows; i++) {
        matrix[rows * i + i] = T(1.0);
    }
}

// Заполняет случайными числами в диапазоне [-1, 1]
template <class T>
void matrixes<T>::MakeRandomNormal() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = static_cast<T>(dis(gen));
    }
}

// Сложение матриц
template <class T>
const matrixes<T> matrixes<T>::operator+(const matrixes<T> &rv) const {
    if ((rows == rv.rows) && (cols == rv.cols)) {
        matrixes<T> mat(rows, cols);

        for (int i = 0; i < rows * cols; i++) {
            mat.matrix[i] = matrix[i] + rv.matrix[i];
        }

        return mat;
    }

    return matrixes<T>(0, 0);
}

// Вычитание матриц
template <class T>
const matrixes<T> matrixes<T>::operator-(const matrixes<T> &rv) const {
    if ((rows == rv.rows) && (cols == rv.cols)) {
        matrixes<T> mat(rows, cols);

        for (int i = 0; i < rows * cols; i++) {
            mat.matrix[i] = matrix[i] - rv.matrix[i];
        }

        return mat;
    }

    return matrixes<T>(0, 0);
}

// Вспомогательное скалярное произведение строки на столбец
template <class T>
const T matrixes<T>::SubstrRowCol(const matrixes<T> &rv, const int n1,
                                    int n2) const {
    T element = T(0.0);

    for (int i = 0; i < cols; i++) {
        element += matrix[cols * n1 + i] * rv.matrix[rv.cols * i + n2];
    }

    return element;
}

// Умножение матриц
template <class T>
const matrixes<T> matrixes<T>::operator*(const matrixes<T> &rv) const {
    T sum = T(0.0);
    int count = -1;

    if (cols != rv.rows) {
        return matrixes<T>(0, 0);
    }

    if (cols == 1) {
        for (int i = 0; i < rows; i++) {
            sum += matrix[i] * rv.matrix[i];
        }

        T *loc = new T[1];
        *loc = sum;
        matrixes<T> mat(loc, 1, 1);

        delete[] loc;
        return mat;
    } else {
        matrixes<T> multmat(rows, rv.cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rv.cols; j++) {
                count++;
                multmat.matrix[count] = SubstrRowCol(rv, i, j);
            }
        }

        return multmat;
    }
}

// Вычисляет определитель через метод Гаусса
template <class T>
T matrixes<T>::GetDet() const {
    if (rows != cols) {
        return static_cast<T>(NAN);
    }

    gauss_method<T> mat(matrix, rows);
    mat.TakeReverse();

    return mat.GetDet();
}
