import matplotlib.pyplot as plt
import math
import Task2 as Task2

#######################################################################
a = -1.0
b = 2.0
N = 20
h = (b - a)/ N
M = 6


def func(x):
    return x**3

def xi(i):
    return a + i*h

x_values = []
y_values = []

for i in xrange(N + 1):
    x_values.append(xi(i))
    y_values.append(func(xi(i)))


def phi_i(i, x):
    return x ** i

def multiply_f(j):
    res = 0
    for i in xrange(0, N + 1):
        res += y_values[i] * phi_i(j, x_values[i])
    return res

def multiply_phi(k, j):
    res = 0
    for i in xrange(0, N + 1):
        res += phi_i(k, x_values[i]) * phi_i(j, x_values[i])
    return res

class Matrixes:
    def __init__(self):
        self.Phi = []
        self.B = []
        for k in xrange(0, N + 1):
            self.Phi.append([])
            for j in xrange(0, N + 1):
                self.Phi[k].append(multiply_phi(k, j))
        for j in xrange(0 , N + 1):
            self.B.append(multiply_f(j))

    def get_matrixes(self):
        return self.Phi, self.B


def multiply(matrix1, matrix2):
        return [[sum(a*b for a,b in zip(matrix1_row,matrix2_col)) for matrix2_col in zip(*matrix2)]
                for matrix1_row in matrix1]


def multiply_vector(matrix, vector):
    if(len(matrix) == len(vector)):
        return [sum(a*b for a,b in zip(matrix_row, vector)) for matrix_row in matrix]

identity_matrix =  lambda n: [[1 if j == i else 0 for j in xrange(n)] for i in xrange(n)]


def get_t_matrix(A, k):
    t = identity_matrix(len(A))
    for p in xrange(len(A)):
        t[p][k] = float(-A[p][k] / A[k][k])
    t[k][k] = (1 / A[k][k])
    return t


def Jourdan_Dauss(A, B):
    for k in xrange(len(A)):
        T = get_t_matrix(A, k)
        A = multiply(T, A)
        B = multiply_vector(T, B)
    return B

matrixes = Matrixes()
A, B = matrixes.get_matrixes()
A_values = []
A_values = Jourdan_Dauss(A, B)

def method_of_least_squares(x):
    res = 0
    for i in xrange(M + 1):
        res += A_values[i] * phi_i(i, x)
    return res

y_least_squares = []
for i in xrange(N + 1):
    y_least_squares.append(method_of_least_squares(x_values[i]))

plt.plot(x_values, y_least_squares)
plt.plot(x_values, y_values)
plt.plot(Task2.x_values, Task2.y_newton)
plt.plot(Task2.x_values, Task2.y_gauss)

plt.show()

##########################norm in L2([a,b]) space###########################

def norm(f):
    res = 0
    for i in xrange(N):
       res += (f(x_values[i + 1]) ** 2 + f(x_values[i]) ** 2) * h / 2
    return math.sqrt(res)

print(abs(norm(func) - norm(Task2.newton)))
print(abs(norm(func) - norm(Task2.gauss)))
print(abs(norm(func) - norm(method_of_least_squares)))
