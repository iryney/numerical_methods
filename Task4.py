from matplotlib import pyplot as plt
import math

a = -1.0
b = 2.0
N = 20
h = (b - a) / N
epsilon = 0.0000001


def func(x):
    return x**3

x_values = []
y_values = []


def xi(i):
    return a + i*h

for i in xrange(0, N+1):
    x_values.append(xi(i))
    y_values.append(func(xi(i)))


def rectangles_method(epsilon):
    iterations = 0
    count = 0
    H = (b - a) / N
    n = N
    res = 0
    res1 = 0
    for i in xrange(0, n):
        res1 += H * func(a + (i + 0.5) * H)
        iterations += 1
    while (math.fabs(res1 - res) / (pow(2,2) - 1)) > epsilon:
        res1 = res
        res = 0
        n *= 2
        H = (b - a) / n
        for i in xrange(0, n):
            res += H * func(a + (i + 0.5) * H)
            iterations += 1
        count += 1
    return res, count, iterations


def trapezoids_method(epsilon):
    iterations = 0
    count = 0
    H = (b - a) / N
    n = N
    res = 0
    res1 = 0
    res1 += (H / 2) * (func(a) + func(b))
    for i in xrange(1, n):
        res1 += H * func(a + H * i)
        iterations += 1
    while (math.fabs(res1 - res) / (pow(2, 2) - 1)) > epsilon:
        res1 = res
        res = 0
        n *= 2
        H = (b - a) / n
        count += 1
        res += (H / 2) * (func(a) + func(b))
        for i in xrange(1, n):
            res += H * func(a + H * i)
            iterations += 1

    return res, count, iterations


def simpsons_method(eps):
    iterations = 0
    count = 0
    H = h
    n = N
    res = 0
    res1 = 0
    for i in xrange(1, n):
        res1 += 2 * func(a + i * H)
        iterations += 1
    for j in xrange(1, n + 1):
        res1 += 4 * func(a + (j - 0.5) * H)
        iterations += 1
    res1 *= H / 6
    while (math.fabs(res1 - res) / (pow(2, 4) - 1)) > eps:
        res1 = res
        res = 0
        n *= 2
        H = (b - a) / n
        count += 1
        res += func(a) + func(b)
        for i in xrange(1, n):
            res += 2 * func(a + i * H)
            iterations += 1
        for j in xrange(1, n + 1):
            res += 4 * func(a + (j - 0.5) * H)
            iterations += 1
        res *= H / 6
    return res, count, iterations


def gauss_method4points(epsilon):
    iterations = 0
    res = 0
    t = [0.861136, 0.339981, -0.339981, -0.861136]
    C = [0.347855, 0.652145, 0.652145, 0.347855]
    for i in xrange(0, 4):
        res += C[i] * func((b + a) / 2 + (b - a) / 2 * t[i])
        iterations += 1
    res *= (b - a) / 2
    return res, iterations


def gauss_method5points(epsilon):
    iterations = 0
    res = 0
    t = [0.90618, 0.538469, 0, -0.538469, -0.90618]
    C = [0.23693, 0.47863, 0.56889, 0.47863, 0.23693]
    for i in xrange(0, 5):
        res += C[i] * func((b + a) / 2 + (b - a) / 2 * t[i])
        iterations += 1
    res *= (b - a) / 2
    return res, iterations


print("Rectangles", rectangles_method(epsilon))
print("Trapezeoids", trapezoids_method(epsilon))
print("Simpsons", simpsons_method(epsilon))
print("Gauss4Points", gauss_method4points(epsilon))
print("Gauss5Points", gauss_method5points(epsilon))

plt.plot(x_values, y_values)
plt.show()


