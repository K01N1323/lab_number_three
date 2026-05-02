#include "../Tests/LinearAlgebraTests.h"
#include "../UI/MatrixUI.h"
#include <iostream>
#include <limits>

using namespace std;

int main() {
    int run_tests = 0;
    while (true) {
        cout << "Запустить тесты перед стартом? (1 - Да, 0 - Нет): ";
        cin >> run_tests;
        if (cin.eof()) {
            return 0;
        }
        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Ошибка: введите 0 или 1\n";
            continue;
        }
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        if (run_tests == 0 || run_tests == 1) {
            break;
        }
        cout << "Ошибка: введите 0 или 1\n";
    }
    if (run_tests == 1) {
        RunAllTests();
    }
    RunUI();
    return 0;
}
