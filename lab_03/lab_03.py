from math import ceil

def dif_matrix_create(table, n):
    for i in range(n):
        tmp = []

        for j in range(n - i):
            tmp.append((table[i + 1][j] - table[i + 1][j + 1]) / (table[0][j] - table[0][i + j + 1]))
        table.append(tmp)

    return table

def choose_points(table, n, x):
    tab_len = len(table[0])
    i_near = min(range(tab_len), key = lambda i: abs(table[0][i] - x)) 
    n_required = ceil(n / 2) 
    
    if (i_near + n_required + 1 > tab_len): 
        i_end = tab_len
        i_start = tab_len - n
    elif (i_near < n_required): 
        i_start = 0
        i_end = n
    else:
        i_start = i_near - n_required + 1
        i_end = i_start + n        

    return [table[0][i_start:i_end], table[1][i_start:i_end]]
    
def newton_interpolation(table, n, x):
    table = choose_points(table, n + 1, x)

    dif_matrix = dif_matrix_create(table, n)
    tmp = 1
    res = 0

    for i in range(n + 1):
        res += tmp * dif_matrix[i + 1][0]
        tmp *= (x - dif_matrix[0][i])

    return res

def f(x):
    return x * x

def xy_table_create(x_beg, step, n):
    x_arr = [(x_beg + i * step) for i in range(n)]
    y_arr = [f(elem) for elem in x_arr]

    return x_arr, y_arr

def xy_table_print(x_arr, y_arr):
    l = len(x_arr)

    for i in range(l):
        print("%.6f %.6f" % (x_arr[i], y_arr[i]))

    print()

def spline_interpolation(x_arr, y_arr, x_value):
    n = len(x_arr)
    i_near = min(range(n), key = lambda i: abs(x_arr[i] - x_value))

    h = [0 if not i else x_arr[i] - x_arr[i - 1] for i in range(n)] # шаг
    
    a_arr = [0 if i < 2 else h[i-1] for i in range(n)]
    b_arr = [0 if i < 2 else -2 * (h[i - 1] + h[i]) for i in range(n)]
    d_arr = [0 if i < 2 else h[i] for i in range(n)]
    f_arr = [0 if i < 2 else -3 * ((y_arr[i] - y_arr[i - 1]) / h[i] - (y_arr[i - 1] - y_arr[i - 2]) / h[i - 1]) for i in range(n)]

    # прямой ход
    ksi = [0 for i in range(n + 1)]
    eta = [0 for i in range(n + 1)]
    for i in range(2, n):
        ksi[i + 1] = d_arr[i] / (b_arr[i] - a_arr[i] * ksi[i])
        eta[i + 1] = (a_arr[i] * eta[i] + f_arr[i]) / (b_arr[i] - a_arr[i] * ksi[i])

    # обратный ход
    c = [0 for i in range(n + 1)]
    for i in range(n - 2, -1, -1):
        c[i] = ksi[i + 1] * c[i + 1] + eta[i + 1]


    a = [0 if i < 1 else y_arr[i-1] for i in range(n)]
    b = [0 if i < 1 else (y_arr[i] - y_arr[i - 1]) / h[i] - h[i] / 3 * (c[i + 1] + 2 * c[i]) for i in range(n)]
    d = [0 if i < 1 else (c[i + 1] - c[i]) / (3 * h[i]) for i in range(n)]


    print(h, '\n', a_arr, '\n', b_arr, '\n', d_arr, '\n', f_arr)
    print()
    print(ksi, '\n', eta)
    print()
    print(a, '\n', b, '\n', c, '\n', d)

    res = a[i_near] + b[i_near] * (x_value - x_arr[i_near - 1]) + c[i_near] * ((x_value - x_arr[i_near - 1]) ** 2) + d[i_near] * ((x_value - x_arr[i_near - 1]) ** 3)

    return res
     
x_beg = float(input("Input beginning value of x: "))
x_step = float(input("Input step for x value: "))
n = int(input("Input points amount: "))

x_tab, y_tab = xy_table_create(x_beg, x_step, n)
print("\nCreated table:")
xy_table_print(x_tab, y_tab)

x = float(input("Input x: "))

res_spline = spline_interpolation(x_tab, y_tab, x)
res_newton = newton_interpolation([x_tab, y_tab], 3, x)
print("\nSpline interpolation: ", res_spline)
print("Newton interpolation: ", res_newton)
print("\nf(x)                : ", f(x))
print("Spline error        : ", abs(f(x) - res_spline))
print("Newton error        : ", abs(f(x) - res_newton), "\n")
