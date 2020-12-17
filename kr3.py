import sys

import numpy as np
import convcode as cc
from sympy import *
import urllib.parse

# 75
MEMORY = np.array([2])
G = np.array([[4, 6, 6]])  # [D^2, D^2+D, D^2+D]
N = '5D^6I^2+D^5I'
M = '5D^5I+3D^3I'
G1 = np.array([[1, 0, 1], [0, 1, 1]])
G2 = np.array([[1, 0, 0, 0], [1, 0, 1, 0], [1, 1, 0, 0]])

# 67
# MEMORY = np.array([2])
# G = np.array([[2, 7, 7]])  # [D, D^2+D+1, D^2+D+1]
# N = 'D^7I^2+D^5I'
# M = '3D^5I+3D^2I'
# G1 = np.array([[1, 1, 0], [1, 1, 1]])
# G2 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]])

WOLFRAM_ALPHA_URL = 'https://www.wolframalpha.com/input/?i='


def format_power(power, negate=False) -> str:
    sign = "-" if negate else ""
    return sign + ("1" if power == 0 else ("D" if power == 1 else "D^" + str(power)))


def format_cell(cell) -> str:
    return "0" if len(cell) == 0 else "".join(cell)


def print_formatted(A, B):
    print("A:")
    for row in A:
        for cell in row:
            print(format_cell(cell), end=" ")
        print()
    print("B:")
    for cell in B:
        print(format_cell(cell), end=" ")
    print()


def get_wolfram_link(input: str, I_replacer='I') -> str:
    return WOLFRAM_ALPHA_URL + urllib.parse.quote_plus(input.replace('I', I_replacer))


if __name__ == '__main__':
    trellis = cc.Trellis(MEMORY, G)
    n = trellis.number_states

    A = [[[] for j in range(n)] for i in range(n)]
    B = [[] for i in range(n)]
    for i in range(n):
        A[i][i].append("1")

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
            power = bin(weight).count('1')
            equation += x * d ** power

            if state_from != 0:
                A[state_to][state_from].append(format_power(power, True))
            else:
                B[state_to].append(format_power(power))
        equations.append(equation)

    print_formatted(A, B)
    g0_function = linsolve(equations, g_all).args[0][0]
    g0_function_wolfram = str(g0_function).replace('**', '^')
    print('dfree:', 'n(D):', 'm(D):', get_wolfram_link('series ' + g0_function_wolfram, '1'), sep='\n')

    print('dC', get_wolfram_link('series (%s) / (1 - (%s))' % (N, M), '*1'), sep='\n')
    # Заменил I на Y, чтобы вольфрам не считал это мнимой единицей
    print('dI', get_wolfram_link('series D[(%s) / (1 - (%s)), Y]' % (N, M), 'Y'), sep='\n')

    print('G:')
    np.savetxt(sys.stdout.buffer, np.kron(G1, G2), fmt='%s', delimiter='')
