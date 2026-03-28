#ifndef MATRIXUI_H
#define MATRIXUI_H

// Выводит меню действий на экран
void print_menu(void);

// Обрабатывает действие, выбранное пользователем
int actions(int flag);

// Запускает бесконечный цикл интерфейса
void run_ui(void);

#endif // MATRIXUI_H
