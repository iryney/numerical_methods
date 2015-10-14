import matplotlib.pyplot as plt
import math
import time


# containers
x_values = []
y_values = []
differences = []

# initial data
a = -1.0
b = 2.0
N = 20
h = (b - a)/ N

def func(x):
    return x ** 3

# find interpolation nodes
def xi(i):
    return a + i*h


# append x_values and y_values
for i in xrange(0, N + 1):
    x_values.append(xi(i))
    y_values.append(func(xi(i)))



def finite_differences(order, i):
   diff = 0.0
   for j in xrange(order + 1):
       diff += ((-1) ** (order - j)) * (math.factorial(order) / (math.factorial(j) *
       math.factorial(order - j)) * func(x_values[i + j]))
   return diff

# get all finite differences
differ = []
differ.append(y_values)
for i in xrange(1, len(y_values)):
    differ.append([])
    for j in xrange(0, len(differ[i - 1]) - 1):
        differ[i].append(differ[i-1][j+1] - differ[i-1][j])


# forward newton
def newton(x):
    t = (x - a) / h
    t_numerator = 1
    Ln = func(a)
    for i in xrange(1, N):
        t_numerator *= t
        Ln += finite_differences(i, 0) * t_numerator / math.factorial(i)
        t -= 1
    return Ln

# first gauss
def gauss(x):
    n = N//2
    x0 = x_values[n]
    t = (x - x0)/ h
    t_numerator = 1
    L = func(x0)
    for i in xrange(1, n+1):
        t_numerator *= (t + i - 1)
        L += t_numerator * finite_differences(2 * i - 1, n - i + 1) / math.factorial(2 * i - 1)
        t_numerator *= (t - i)
        L += t_numerator * finite_differences(2 * i, n - i) / math.factorial(2 * i)
    return L


y_newton = []
y_gauss = []


start_time_1 = time.clock()
for i in xrange(0, N + 1):
    y_newton.append(newton(xi(i)))
    y_gauss.append(gauss(xi(i)))
print (time.clock() - start_time_1)


plt.plot(x_values, y_values)
plt.plot(x_values, y_newton)
plt.plot(x_values, y_gauss)

plt.show()


