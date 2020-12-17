import numpy as np
import convcode as cc
from sympy import *
import urllib.parse


MEMORY = np.array([2])
G = np.array([[4, 6, 6]])  # [D^2, D^2+D, D^2+D]
G1 = np.array([[1, 1, 0], [1, 1, 1]])
G2 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]])

WOLFRAM_ALPHA_URL = 'https://www.wolframalpha.com/input/?i='


def get_wolfram_link(input: str) -> str:
    return WOLFRAM_ALPHA_URL + urllib.parse.quote_plus(input)


if __name__ == '__main__':
    trellis = cc.Trellis(MEMORY, G)
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

    g0_function = linsolve(equations, g_all).args[0][0]
    # Заменил I на Y, чтобы вольфрам не считал это мнимой единицей
    g0_function_wolfram = str(g0_function).replace('**', '^').replace('I', 'Y')

    print(equations)
    print(get_wolfram_link("series " + g0_function_wolfram))
    # print(get_wolfram_link("series D[%s, Y]" % g0_function_wolfram))
    print(np.kron(G1, G2))

# D^2 + D
# 0 1 1 -> 3
