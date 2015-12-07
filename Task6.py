from matplotlib import pyplot as plt
import math

a = 0
b = 2 * math.pi
# a = 0.7
# b = 1.0
alpha0 = 1
# alpha0 = 0
alpha1 = 0
# alpha1 = 1
betha0 = 1
# betha0 = 1
betha1 = 0
# betha1 = 0
A = 0
B = 0
# A = 1
# B = 1

def f(x):
    return math.cos(x)
    #return 1.5
def p(x):
    # return 0
    return 1
def q(x):
    return 1

N = 40
h = (b - a) / N
x_values = []
x_values.append(a)
for i in xrange(1, N + 1):
    x_values.append(a + h * i)


M = [0.0 for i in range(N)]
for i in xrange(1, N):
    mi = ((h ** 2) * q(x_values[i]) - 2) / (1 + (h / 2) * p(x_values[i]))
    M[i] = mi

K = [0.0 for i in range(N)]
for i in xrange(1, N):
    ki = (1 - (h / 2) * p(x_values[i])) / (1 + (h / 2) * p(x_values[i]))
    K[i] = ki

F = [0.0 for i in range(N)]
for i in xrange(1, N):
    fi = f(x_values[i]) / (1 + (h / 2) * p(x_values[i]))
    F[i] = fi


C = []
C.append((alpha1 / (h * alpha0 - alpha1)))
for i in xrange(1, N):
    ci = 1 / (M[i] - K[i] * C[i - 1])
    C.append(ci)

D = []
D.append(A * h / (h * alpha0 - alpha1))
for i in xrange(1, N):
    di = C[i] * ((h ** 2) * F[i] - K[i] * D[i - 1])
    D.append(di)

y_values = [0.0 for i in xrange(N + 1)]
Yn = (B * h + betha1 * D[N - 1]) / (betha0 * h + betha1 * (C[N - 1] + 1))
y_values[N] = Yn
for i in xrange(N - 1, 0, -1):
    y_values[i] = D[i] - C[i] * y_values[i + 1]



result = zip(x_values, y_values)
count = 0
for item in result:
    print count, item
    count += 1


plt.plot(x_values, y_values)
plt.show()