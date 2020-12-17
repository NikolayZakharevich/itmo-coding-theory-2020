import numpy as np
import convcode as cc
from sympy import *

# memory = np.array([2])
# g_matrix = np.array([[1, 3, 7]])

# memory = np.array([2])
# g_matrix = np.array([[5, 7]])

memory = np.array([3])
g_matrix = np.array([[13, 11]])

trellis = cc.Trellis(memory, g_matrix)
n = trellis.number_states

trellis.visualize()

equations = []
g_all = symbols('g0:%d' % n)
d = symbols('D')
i = symbols('I')
for state_to in range(n):
    equation = -g_all[state_to]
    for state_from in range(n):
        if not np.count_nonzero(trellis.next_state_table[state_from, :] == state_to) or (
                state_to == 0 and state_from == 0):
            continue
        bit = np.where(trellis.next_state_table[state_from, :] == state_to)[0][0]
        weight = int(trellis.output_table[state_from, bit])
        x = 1
        if state_from != 0:
            x = g_all[state_from]
        if bit == 1:
            x *= i
        equation += x * d ** bin(weight).count('1')

    equations.append(equation)

print(equations)
g0_function = linsolve(equations, g_all).args[0][0]

print(g0_function)
