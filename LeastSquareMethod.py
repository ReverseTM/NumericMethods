import sys
import matplotlib.pyplot as plt
from GaussMethod import *


# Функция для нахождения коэффициентов многочлена МНК
def least_squares_polynomial(x, y, degree):
    a = np.zeros((degree + 1, degree + 1))
    b = np.zeros(degree + 1)

    # Заполняем матрицу A и вектор B
    for i in range(degree + 1):
        for j in range(degree + 1):
            a[i, j] = np.sum(x**(i + j))
        b[i] = np.sum(y * x**i)

    # Решаем нормальную систему уравнений для нахождения коэффициентов
    coefficients = solve_linear_system(a, b)
    return coefficients


# Функция для вычисления значения многочлена
def evaluate_polynomial(coefficients, x):
    degree = len(coefficients) - 1
    result = np.zeros_like(x)

    # Считаем значение многочлена
    for i in range(degree + 1):
        result += coefficients[i] * (x**i)

    return result


# Функция для вычисления суммы квадратов ошибок
def sum_of_squares_error(y, y_approx):
    return np.sum((y - y_approx)**2)


def main():
    # Таблично заданные данные
    x_i = np.array([0.1, 0.5, 0.9, 1.3, 1.7, 2.1])
    y_i = np.array([-2.2026, -0.19315, 0.79464, 1.5624, 2.2306, 2.8419])

    # Для первой степени
    degree_1 = 1
    coefficients_1 = least_squares_polynomial(x_i, y_i, degree_1)
    y_approx_1 = evaluate_polynomial(coefficients_1, x_i)
    error_1 = sum_of_squares_error(y_i, y_approx_1)

    # Для второй степени
    degree_2 = 2
    coefficients_2 = least_squares_polynomial(x_i, y_i, degree_2)
    y_approx_2 = evaluate_polynomial(coefficients_2, x_i)
    error_2 = sum_of_squares_error(y_i, y_approx_2)

    # Открываем файл для записи
    with open('result.txt', 'w') as file:
        sys.stdout = file

        # Вывод результатов
        print(f'Coefficients of a polynomial of the 1st degree: {coefficients_1}')
        print(f'Sum of squared errors for 1st degree: {error_1}')
        print(f'Coefficients of a polynomial of the 2nd degree: {coefficients_2}')
        print(f'Sum of squared errors for 2nd degree: {error_2}')

    # Построение графиков
    x_values = np.linspace(min(x_i), max(x_i), 100)
    y_true = np.interp(x_values, x_i, y_i)

    plt.scatter(x_i, y_i, label='Табличные данные')
    plt.plot(x_values, evaluate_polynomial(coefficients_1, x_values), label=f'Полином 1-ой степени')
    plt.plot(x_values, evaluate_polynomial(coefficients_2, x_values), label=f'Полином 2-ой степени')
    plt.plot(x_values, y_true, '--', label='Приближаемая функция')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
    