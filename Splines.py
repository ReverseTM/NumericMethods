import matplotlib.pyplot as plt
import numpy as np


def solution(matrix: [[float]], vector_b: [float]) -> [float]:
    vector_alphas = [-matrix[0][1] / matrix[0][0]]
    vector_betas = [vector_b[0] / matrix[0][0]]

    for i in range(1, len(matrix) - 1):
        y_i = matrix[i][i] + matrix[i][i - 1] * vector_alphas[i - 1]
        a_i = -matrix[i][i + 1] / y_i
        b_i = (vector_b[i] - matrix[i][i - 1] * vector_betas[i - 1]) / y_i
        vector_alphas.append(a_i)
        vector_betas.append(b_i)

    y_n = (matrix[-1][-1] + matrix[-1][-2] * vector_alphas[-1])
    vector_x = [round((vector_b[-1] - matrix[-1][-2] * vector_betas[-1]) / y_n, 6)]

    for i in range(len(matrix) - 2, -1, -1):
        vector_x.insert(0, round(vector_alphas[i] * vector_x[0] + vector_betas[i], 6))

    return vector_x


def get_index(vector, x):
    for i in range(len(vector) - 1):
        if vector[i] <= x <= vector[i + 1]:
            return i
    return 0


def get_spline_matrix(x_i: [float], f_i: [float]) -> [[float]]:
    res = [[0 for i in range(len(x_i) - 2)] for j in range(len(x_i) - 2)]
    vector = [0] * (len(f_i) - 2)
    h_i = [x_i[i] - x_i[i - 1] for i in range(1, len(x_i))]
    res[0][0] = 2 * (h_i[0] + h_i[1])
    res[0][1] = h_i[1]
    vector[0] = 3 * ((f_i[2] - f_i[1]) / h_i[1] - (f_i[1] - f_i[0]) / h_i[0])
    j = 1
    for i in range(2, len(x_i) - 2):
        res[i - 1][i - 2] = h_i[i]
        res[i - 1][i - 1] = 2 * (h_i[i] + h_i[i + 1])
        res[i - 1][i] = h_i[i]
        vector[j] = 3 * ((f_i[i + 1] - f_i[i]) / h_i[i] - (f_i[i] - f_i[i - 1]) / h_i[i - 1])
        j += 1
    res[2][1] = h_i[2]
    res[2][2] = 2 * (h_i[2] + h_i[3])
    vector[2] = 3 * ((f_i[4] - f_i[3]) / h_i[3] - (f_i[3] - f_i[2]) / h_i[2])
    return res, vector


def get_spline(x_i, f_i, x):
    spline_matrix_and_vector = get_spline_matrix(x_i, f_i)
    h_i = [x_i[i] - x_i[i - 1] for i in range(1, len(x_i))]
    c_i = [0] + solution(*spline_matrix_and_vector)
    a_i = [f_i[i] for i in range(0, len(x_i) - 1)]
    d_i, b_i = [], []
    for i in range(1, len(x_i) - 1):
        b_i.append((f_i[i] - f_i[i - 1]) / h_i[i - 1] - 1 / 3 * (c_i[i] + 2 * c_i[i - 1]) * h_i[i - 1])
        d_i.append((c_i[i] - c_i[i - 1]) / (3 * h_i[i - 1]))
    b_i.append((f_i[4] - f_i[3]) / h_i[3] - 2 / 3 * h_i[3] * c_i[3])
    d_i.append(-c_i[3] / (3 * h_i[3]))
    index = get_index(x_i, x)
    res = a_i[index] + b_i[index] * (x - x_i[index]) + c_i[index] * (x - x_i[index]) ** 2 + d_i[
        index] * (x - x_i[index]) ** 3
    return res, a_i, b_i, c_i, d_i


def main():
    x_i = [0.1, 0.5, 0.9, 1.3, 1.7]
    f_i = [-2.2026, -0.19315, 0.79464, 1.5624, 2.2306]

    x0 = 0.89

    res, a, b, c, d = get_spline(x_i, f_i, x0)
    print(f'coefficient a = {a}')
    print(f'coefficient b = {b}')
    print(f'coefficient c = {c}')
    print(f'coefficient d = {d}')
    print(f"Spline value in the area {x0} = {res}")

    x_points = np.linspace(x_i[0], x_i[-1], 1000)

    # Генерируем значения сплайна для построения графика
    # y_spline = [get_spline(x_i, f_i, x)[0] for x in x_points]
    # y_spline1 = [get_spline(x_i, f_i, x)[0] for x in x_i]
    # 
    # # Строим график
    # plt.plot(x_points, y_spline, label='Cubic Spline')
    # plt.scatter(x_i, y_spline1, c='red', label='Spline Points')
    # plt.legend()
    # plt.title('Cubic Spline Function')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.show()


if __name__ == '__main__':
    main()
    