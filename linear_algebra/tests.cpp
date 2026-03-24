#include <chrono>
#include <iostream>
#include <random>

#include "LU_decompozition.h"
#include "default_gauss.h"
#include "gauss_method.h"
#include "matrixes.h"
#include "tests.h"

using namespace std;

struct ExecutionTime {
    double lu_time;
    double solve_time;
};

// Замеряет время решения случайной СЛАУ методом Гаусса без ведущего элемента
double Solve_with_default_gauss(int n) {

    matrixes<double> Matrix(n, n);
    Matrix.MakeRandomNormal();

    double *b = new double[n];

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    for (int i = 0; i < n; i++) {
        b[i] = dis(gen);
    }

    default_gauss<double> Matrix_for_gauss(Matrix.GetMatrix(), n);

    auto start_gauss = std::chrono::high_resolution_clock::now();
    Matrix_for_gauss.SolveForTests(b);
    auto end_gauss = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time = end_gauss - start_gauss;

    delete[] b;
    return time.count();
}

// Замеряет время решения СЛАУ методом Гаусса с ведущим элементом
double Solve_with_gauss(int n) {

    matrixes<double> Matrix(n, n);
    Matrix.MakeRandomNormal();

    double *b = new double[n];

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    for (int i = 0; i < n; i++) {
        b[i] = dis(gen);
    }

    gauss_method<double> Matrix_for_gauss(Matrix.GetMatrix(), n);

    auto start_gauss = std::chrono::high_resolution_clock::now();
    Matrix_for_gauss.SolveForTests(b);
    auto end_gauss = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time = end_gauss - start_gauss;

    delete[] b;
    return time.count();
}

// Раздельно замеряет время факторизации и решения для метода LU-разложения
ExecutionTime Solve_with_LU(int n) {
    matrixes<double> Matrix(n, n);
    Matrix.MakeRandomNormal();

    double *b = new double[n];

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    for (int i = 0; i < n; i++) {
        b[i] = dis(gen);
    }

    auto start_lu = std::chrono::high_resolution_clock::now();
    LU_Decomposition<double> Matrix_for_LU(Matrix.GetMatrix(), n);
    auto end_lu = std::chrono::high_resolution_clock::now();

    auto start_solve = std::chrono::high_resolution_clock::now();
    Matrix_for_LU.SolveForTests(b);
    auto end_solve = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff_lu = end_lu - start_lu;
    std::chrono::duration<double> diff_solve = end_solve - start_solve;

    delete[] b;
    return {diff_lu.count(), diff_solve.count()};
}

// Тестирование скорости решения k систем с одинаковой матрицей (LU против
// Гаусса)
ExecutionTime Solving_immutable_b(int k) {
    double summ_time_for_gauss = 0;
    double summ_time_for_lu = 0;
    int n = 500;

    matrixes<double> matrix(n, n);
    matrix.MakeRandomNormal();

    gauss_method<double> matrix_for_gauss(matrix.GetMatrix(), n);

    auto start_lu_decompoze = std::chrono::high_resolution_clock::now();
    LU_Decomposition<double> matrix_for_LU(matrix.GetMatrix(), n);
    auto end_lu_decompoze = std::chrono::high_resolution_clock::now();

    summ_time_for_lu += (end_lu_decompoze - start_lu_decompoze).count();

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    double *b = new double[n];

    for (int i = 0; i < k; i++) {
        for (int k = 0; k < 500; k++) {
            b[k] = dis(gen);
        }

        auto start_gauss = std::chrono::high_resolution_clock::now();
        matrix_for_gauss.SolveForTests(b);
        auto end_gauss = std::chrono::high_resolution_clock::now();

        auto start_lu = std::chrono::high_resolution_clock::now();
        matrix_for_LU.SolveForTests(b);
        auto end_lu = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> time_for_gauss = end_gauss - start_gauss;
        std::chrono::duration<double> time_for_lu = end_lu - start_lu;

        summ_time_for_gauss += time_for_gauss.count();
        summ_time_for_lu += time_for_lu.count();
    }

    delete[] b;

    return {summ_time_for_gauss, summ_time_for_lu};
}

// Тестирование погрешностей для матрицы Гильберта обычным Гауссом
ExecutionTime Gilbert_tests_for_default_gauss(int n) {
    matrixes<double> matrix(n, n);
    matrix.MakeGilbert();

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    double *b = new double[n];

    for (int i = 0; i < n; i++) {
        b[i] = 0;

        for (int j = 0; j < n; j++) {
            b[i] += matrix.GetIJ(i, j);
        }
    }

    default_gauss<double> matrix_for_gauss(matrix.GetMatrix(), n);
    double *x = matrix_for_gauss.Solve(b);

    double result_for_normal = 0;

    for (int i = 0; i < n; i++) {
        result_for_normal += std::abs(x[i] - 1.0);
    }

    double result_for_discrepancy = 0.;

    for (int i = 0; i < n; i++) {
        double Hx_i = 0;

        for (int j = 0; j < n; j++) {
            Hx_i += matrix.GetIJ(i, j) * x[j];
        }

        result_for_discrepancy += std::abs(b[i] - Hx_i);
    }

    delete[] b;

    return {result_for_normal, result_for_discrepancy};
}

// Тестирование погрешностей для матрицы Гильберта Гауссом с ведущим элементом
ExecutionTime Gilbert_tests_for_gauss(int n) {
    matrixes<double> matrix(n, n);
    matrix.MakeGilbert();

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    double *b = new double[n];

    for (int i = 0; i < n; i++) {
        b[i] = 0;

        for (int j = 0; j < n; j++) {
            b[i] += matrix.GetIJ(i, j);
        }
    }

    gauss_method<double> matrix_for_gauss(matrix.GetMatrix(), n);
    double *x = matrix_for_gauss.Solve(b);

    double result_for_normal = 0;

    for (int i = 0; i < n; i++) {
        result_for_normal += std::abs(x[i] - 1.0);
    }

    double result_for_discrepancy = 0.;

    for (int i = 0; i < n; i++) {
        double Hx_i = 0;

        for (int j = 0; j < n; j++) {
            Hx_i += matrix.GetIJ(i, j) * x[j];
        }

        result_for_discrepancy += std::abs(b[i] - Hx_i);
    }

    delete[] b;

    return {result_for_normal, result_for_discrepancy};
}

// Тестирование погрешностей для матрицы Гильберта LU-разложением
ExecutionTime Gilbert_tests_for_lu(int n) {
    matrixes<double> matrix(n, n);
    matrix.MakeGilbert();

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    double *b = new double[n];

    for (int i = 0; i < n; i++) {
        b[i] = 0;

        for (int j = 0; j < n; j++) {
            b[i] += matrix.GetIJ(i, j);
        }
    }

    LU_Decomposition<double> matrix_for_lu(matrix.GetMatrix(), n);
    double *x = matrix_for_lu.Solve(b);

    double result_for_normal = 0;

    for (int i = 0; i < n; i++) {
        result_for_normal += std::abs(x[i] - 1.0);
    }

    double result_for_discrepancy = 0.;

    for (int i = 0; i < n; i++) {
        double Hx_i = 0;

        for (int j = 0; j < n; j++) {
            Hx_i += matrix.GetIJ(i, j) * x[j];
        }

        result_for_discrepancy += std::abs(b[i] - Hx_i);
    }

    delete[] b;

    return {result_for_normal, result_for_discrepancy};
}
// Главная функция вызова всех тестов и вывода результатов на консоль
void run_all_tests() {
    cout << "НАЧАЛО ТЕСТОВ МЕТОДОВ РЕШЕНИЯ СЛАУ" << endl;

    cout << "МЕТОД ГАУССА БЕЗ ВЫБОРА ВЕДУЩЕГО ЭЛЕМЕНТА" << endl;
    cout << "Время нахождения решения для матрицы 100 на 100" << endl;
    cout << Solve_with_default_gauss(100) << endl;
    cout << "Время нахождения решения для матрицы 200 на 200" << endl;
    cout << Solve_with_default_gauss(200) << endl;
    cout << "Время нахождения решения для матрицы 500 на 500" << endl;
    cout << Solve_with_default_gauss(500) << endl;

    cout << "МЕТОД ГАУССА С ВЫБОРОМ ВЕДУЩЕГО ЭЛЕМЕНТА" << endl;
    cout << "Время нахождения решения для матрицы 100 на 100" << endl;
    cout << Solve_with_gauss(100) << endl;
    cout << "Время нахождения решения для матрицы 200 на 200" << endl;
    cout << Solve_with_gauss(200) << endl;
    cout << "Время нахождения решения для матрицы 500 на 500" << endl;
    cout << Solve_with_gauss(500) << endl;

    cout << "LU-ДЕКОМПОЗИЦИЯ" << endl;

    cout << "Время нахождения решения для матрицы 100 на 100" << endl;
    ExecutionTime t100 = Solve_with_LU(100);
    cout << "Время LU: " << t100.lu_time << endl;
    cout << "Время Solve: " << t100.solve_time << endl;
    cout << "Суммарно: " << t100.lu_time + t100.solve_time << endl;

    cout << "Время нахождения решения для матрицы 200 на 200" << endl;
    ExecutionTime t200 = Solve_with_LU(200);
    cout << "Время LU: " << t200.lu_time << endl;
    cout << "Время Solve: " << t200.solve_time << endl;
    cout << "Суммарно: " << t200.lu_time + t200.solve_time << endl;

    cout << "Время нахождения решения для матрицы 500 на 500" << endl;
    ExecutionTime t500 = Solve_with_LU(500);
    cout << "Время LU: " << t500.lu_time << endl;
    cout << "Время Solve: " << t500.solve_time << endl;
    cout << "Суммарно: " << t500.lu_time + t500.solve_time << endl;

    cout << "Тесты скорости решения СЛАУ LU разложением при неизменной матрице "
            "коэффициентов"
         << endl;
    cout << "Для k == 1" << std::endl;
    ExecutionTime t3 = Solving_immutable_b(1);
    cout << "Время для метода Гаусса: " << t3.lu_time << endl;
    cout << "Время для LU разложения: " << t3.solve_time << endl;
    cout << endl;

    cout << "Для k == 10" << std::endl;
    ExecutionTime t10 = Solving_immutable_b(10);
    cout << "Время для метода Гаусса: " << t10.lu_time << endl;
    cout << "Время для LU разложения: " << t10.solve_time << endl;
    cout << endl;

    cout << "Для k == 100" << std::endl;
    ExecutionTime t100_ = Solving_immutable_b(100);
    cout << "Время для метода Гаусса: " << t100_.lu_time << endl;
    cout << "Время для LU разложения: " << t100_.solve_time << endl;
    cout << endl;

    cout << "Тестирование невязки и относительной погрешности для Метода Гаусса"
         << endl;
    cout << "Для n = 5" << endl;
    ExecutionTime t5 = Gilbert_tests_for_default_gauss(5);
    cout << "Отностительная погрешность: " << t5.lu_time << endl;
    cout << "Невязка: " << t5.solve_time << endl;

    cout << "Для n = 10" << endl;
    ExecutionTime t10_ = Gilbert_tests_for_default_gauss(10);
    cout << "Отностительная погрешность: " << t10_.lu_time << endl;
    cout << "Невязка: " << t10_.solve_time << endl;

    cout << "Для n = 10" << endl;
    ExecutionTime t15 = Gilbert_tests_for_default_gauss(15);
    cout << "Отностительная погрешность: " << t15.lu_time << endl;
    cout << "Невязка: " << t15.solve_time << endl;

    cout
        << "Тестирование невязки и относительной погрешности для Метода Гаусса "
           "c выбором ведушего элемента"
        << endl;
    cout << "Для n = 5" << endl;
    t5 = Gilbert_tests_for_gauss(5);
    cout << "Отностительная погрешность: " << t5.lu_time << endl;
    cout << "Невязка: " << t5.solve_time << endl;

    cout << "Для n = 10" << endl;
    t10_ = Gilbert_tests_for_gauss(10);
    cout << "Отностительная погрешность: " << t10_.lu_time << endl;
    cout << "Невязка: " << t10_.solve_time << endl;

    cout << "Для n = 10" << endl;
    t15 = Gilbert_tests_for_gauss(15);
    cout << "Отностительная погрешность: " << t15.lu_time << endl;
    cout << "Невязка: " << t15.solve_time << endl;

    cout
        << "Тестирование невязки и относительной погрешности для Lu разложения "
        << endl;
    cout << "Для n = 5" << endl;
    t5 = Gilbert_tests_for_lu(5);
    cout << "Отностительная погрешность: " << t5.lu_time << endl;
    cout << "Невязка: " << t5.solve_time << endl;

    cout << "Для n = 10" << endl;
    t10_ = Gilbert_tests_for_lu(10);
    cout << "Отностительная погрешность: " << t10_.lu_time << endl;
    cout << "Невязка: " << t10_.solve_time << endl;

    cout << "Для n = 10" << endl;
    t15 = Gilbert_tests_for_lu(15);
    cout << "Отностительная погрешность: " << t15.lu_time << endl;
    cout << "Невязка: " << t15.solve_time << endl;
}
