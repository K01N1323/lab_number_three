#include "../Tests/LinearAlgebraTests.h"
#include "../UI/MatrixUI.h"
#include <iostream>

using namespace std;

int main() {
    int run_tests;
    cout << "Запустить тесты перед стартом? (1 - Да, 0 - Нет): ";
    if (cin >> run_tests && run_tests == 1) {
        run_all_tests();
    }
    run_ui();
    return 0;
}
