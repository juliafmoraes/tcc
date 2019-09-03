#!/usr/bin/env python
from methods import TDMA_solver
#
def test_1():
    """
    https://www.ufrgs.br/reamat/CalculoNumerico/livro-sci/sdsl-metodo_da_matriz_tridiagonal.html
    solution is (1,2,-1,0,1)
    :return:
    """
    print("starting test_1")
    diag = [2] * 5
    l_diag = [0, 1, 1, 1, 1]
    u_diag = [1, 1, 1, 1, 0]
    d = [4, 4, 0, 0, 2]

    solution = TDMA_solver(diag, l_diag, u_diag, d)
    print(solution)
    return solution


if __name__=="__main__":
    print("unit tests for tdma-solver")
    x = test_1()