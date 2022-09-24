# Аппроксимация ф-и
# Наилучшее среднеквадратичное значение
import matplotlib.pyplot as plt
import numpy as np

def f(x_arr, coef):
    res = np.zeros(len(x_arr))

    n = len(coef)
    for i in range(n):
        res += coef[i] * (x_arr**i)

    return res                   

def file_read(filename):
    f = open(filename, "r")
    x = []
    y = []
    p = []

    for line in f:
        line = line.split(" ")
        x.append(float(line[0]))
        y.append(float(line[1]))
        p.append(float(line[2]))

    return x, y, p

def table_print(x, y, p):
    n = len(x)
    print("%8s |%8s |%8s " % ("x   ", "y   ", "p   "))

    for i in range(n):
        print("%.6f |%.6f |%.6f" % (x[i], y[i], p[i]))
    print()

def matrix_print(matrix):
    for el in matrix:
        print(el)

def gauss_method(matrix):
    n = len(matrix)
    # приводим к треугольному виду
    for k in range(n):
        for i in range(k + 1, n):
            coef = -(matrix[i][k] / matrix[k][k])
            for j in range(k, n + 1):
                matrix[i][j] += coef * matrix[k][j]

    #print("\ntriangled:")
    #matrix_print(matrix)

    # находим неизвестные
    a = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            matrix[i][n] -= a[j] * matrix[i][j]
        a[i] = matrix[i][n] / matrix[i][i]

    return a

### Вычислить значение
def solve(x, y, p, n): #n - кол-во искомых коэффициентов
    l = len(x)

    sum_x_n = [sum([x[i]**j * p[i] for i in range(l)]) for j in range(n*2 -1)]
    sum_x_n = [round(el, 6) for el in sum_x_n]
    sum_yx_n = [sum([x[i]**j * p[i] * y[i] for i in range(l)]) for j in range(n)]

    matrix = [sum_x_n[i:i+n] for i in range(n)]

    for i in range(n):
        matrix[i].append(round(sum_yx_n[i], 6))
        #matrix[i].append(sum_yx_n[i])
    matrix_print(matrix)

    return gauss_method(matrix)
    

### Отобразить результат
def figure_show(a, x, y, p, n): #, a1):
    colors = ['b', 'g', 'c', 'm', 'y']
    #lbl = ['1', '2', '3', '4', '5']

    t = np.arange(x[0] - 2, x[len(x) - 1] + 2, 0.01)

    plt.figure(1)
    
    plt.xlabel("x")
    plt.ylabel("y")

    for i in range(n):
        plt.plot(t, f(t, a[i]), colors[i], label = 'n = ' + str(i + 1)) #, label = 'различные веса')
    #plt.plot(t, f(t, a[1]), 'g')
    #plt.plot(t, f(t, a[2]), 'y')
    #plt.plot(t, f(t, a1), 'k--', label = 'p_i = 1')

    plt.legend()

    for i in range(len(x)):
        plt.plot(x[i], y[i], 'ro', markersize = p[i] + 2)

    plt.show()

if __name__ == '__main__':
    x, y, p = file_read("data.txt")
    N = len(x) - 1 # Степень многочлена    #2
    n = int(input("Введите порядок полинома: "))
    table_print(x, y, p)

    a = []
    for i in range(1, n + 1):
        a.append(solve(x, y, p, i + 1))

    #p1 = [1 for i in range(len(x))]
    #a1 = solve(x, y, p1, n + 1)
    #print("\na:", a)
    figure_show(a, x, y, p, n)  #, a1)