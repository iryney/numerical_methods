from matplotlib import pyplot as plt
import math

a = 0.20
b = 1.20
y0 = 0.25
epsilon = 0.0001


def func(x, y):
    return 0.188 * (x ** 2 + math.sin(1.5 * x)) + 0.885 * y

def euler_method(N):
    x_values = []
    y_values = []
    h = (b - a) / N
    x_values.append(a)
    y_values.append(y0)
    for i in xrange(1, N):
        x_values.append(a + i * h)
        Yk = y_values[i - 1] + h * func(x_values[i - 1], y_values[i - 1])
        Yk1 = y_values[i - 1] + h * (func(x_values[i - 1], y_values[i - 1]) + func(x_values[i - 1], Yk)) / 2
        while(math.fabs(Yk1 - Yk) > epsilon):
            Yk = y_values[i - 1] + h * (func(x_values[i - 1], y_values[i - 1]) + func(x_values[i - 1], Yk1)) / 2
            Yk1 = y_values[i - 1] + h * (func(x_values[i - 1], y_values[i - 1]) + func(x_values[i - 1], Yk)) / 2
        y_values.append(Yk1)
    return x_values, y_values

def runge_kutta_methods(N):
    x_values = []
    y_values = []
    h = (b - a) / N
    x_values.append(a)
    y_values.append(y0)
    for i in xrange(1, N):
        x_values.append(a + h * i)
        k1 = func(x_values[i - 1], y_values[i - 1])
        k2 = func(x_values[i - 1] + (1 / 2) * h, y_values[i - 1] + (h / 2) * k1)
        k3 = func(x_values[i - 1] + (1 / 2) * h, y_values[i - 1] + (h / 2) * k2)
        k4 = func(x_values[i - 1] + h, y_values[i - 1] + h * k3)
        yi = y_values[i - 1] + (h / 6) * (k1 + (2 * k2) + (2 * k3) + k4)
        y_values.append(yi)
    return x_values, y_values


def gauss_method5points(s, f):
    res = 0
    t = [0.906, 0.538469, 0, -0.538469, -0.90618]
    C = [0.23693, 0.47863, 0.56889, 0.47863, 0.23693]
    for i in xrange(0, 5):
        res += C[i] * f(s, (1.0 + 0.0) / 2 + (1.0 - 0.0) / 2 * t[i])
    res *= (1.0 - 0.0) / 2
    return res


def alpha_extrapolation(s, t):
    res = 1.0
    if s == 0:
        return 1
    for i in xrange(0, s):
       res *= (t + i)
    res /= math.factorial(s)
    return res


def differences(N):
    x_values, y_values = runge_kutta_methods(N)
    differ = []
    func_values = []
    for i in xrange(len(x_values)):
        func_values.append(func(x_values[i], y_values[i]))

    differ.append(func_values)
    for i in xrange(1, len(func_values)):
        differ.append([])
        for j in xrange(0, len(differ[i - 1]) - 1):
            differ[i].append(differ[i-1][j+1] - differ[i-1][j])
    return differ



def extrapolation_method(N, M):
    differ = differences(N)
    x_values = []
    y_values = []
    h = (b - a) / N
    x_values.append(a)
    y_values.append(y0)
    for i in xrange(1, N):
        x_values.append(a + h * i)
        adds = 0
        for s in xrange(0, M + 1):
            if(s <= (i - 1)):
                adds += gauss_method5points(s, alpha_extrapolation) * differ[s][i - 1 - s]
            else:
                break
        yi = y_values[i - 1] + h * adds
        y_values.append(yi)
    return x_values, y_values


def alpha_interpolation(s, t):
    res = 1.0
    if s == 0:
        return 1
    for i in xrange(0, s):
       res *= (t + i - 1)
    res /= math.factorial(s)
    return res

def interpolation_method(N, M):
    differ = differences(N)
    x_values = []
    y_values = []
    h = (b - a) / N
    x_values.append(a)
    y_values.append(y0)
    for i in xrange(1, N):
        x_values.append(a + h * i)
        adds = 0
        for s in xrange(1, M + 1):
            if(s <= i):
                adds += gauss_method5points(s, alpha_interpolation) * differ[s][i - s]
            else:
                break
        delta = y_values[i - 1] + h * adds
        y1 = h * func(x_values[i], y_values[i - 1]) + delta
        y2 = h * func(x_values[i], y1) + delta
        while(math.fabs(y1 - y2) > epsilon):
            y1 = h * func(x_values[i], y2)
            y2 = h * func(x_values[i], y1)
        y_values.append(y1)
    return x_values, y_values


x_euler, y_euler = euler_method(100)
x_runge, y_runge = runge_kutta_methods(100)
x_ex, y_ex = extrapolation_method(100, 49)
x_int, y_int = interpolation_method(100, 49)

plt.plot(x_runge, y_runge)
plt.plot(x_euler, y_euler)
plt.plot(x_ex, y_ex)
plt.plot(x_ex, y_ex)
plt.show()