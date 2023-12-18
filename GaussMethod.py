import numpy as np


# Функция для решения системы линейных уравнений A * x = b методом Гаусса
def solve_linear_system(a, b):
    n = len(a)
    # Создаем копии матрицы A и вектора b
    u = a.copy()
    y = b.copy()

    # Прямой ход метода Гаусса
    for k in range(n):
        # Находим ведущий элемент
        pivot = u[k, k]
        # Обновляем строку k матрицы U и соответствующий элемент вектора y
        u[k, :] /= pivot
        y[k] /= pivot
        for i in range(k + 1, n):
            factor = u[i, k]
            u[i, :] -= factor * u[k, :]
            y[i] -= factor * y[k]

    # Обратный ход метода Гаусса
    x = np.zeros(n)
    for k in range(n - 1, -1, -1):
        x[k] = y[k]
        for i in range(k + 1, n):
            x[k] -= u[k, i] * x[i]

    return x
