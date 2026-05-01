#include "../Tests/LinearAlgebraTests.h"
#include "../UI/MatrixUI.h"
#include <iostream>
#include <limits>

using namespace std;

int main() {
    int run_tests;
    cout << "Запустить тесты перед стартом? (1 - Да, 0 - Нет): ";
    if (!(cin >> run_tests)) {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        run_tests = 0;
    } else {
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    if (run_tests == 1) {
        RunAllTests();
    }
    RunUI();
    return 0;
}
