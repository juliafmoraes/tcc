#!/usr/bin/env python
import numpy as np
import pandas as pd
import openpyxl

from simulations.simulation import Simulation
from methods import TDMA_solver


class ClassicFick(Simulation):

    def __init__(self, dt, dx, T, n, D: float, lower, upper):
        """
        :type D: float
        """
        super().__init__(dt, dx, T, n)
        self.D_coef = D
        self.boundary_conditions(lower, upper)

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper


def create_d(c, a, lb):
    d = [-1 * element for element in c]
    d[0] = d[0] - a * lb
    return d


def run():
    # condicoes de teste
    lower_bound = 5
    upper_bound = 0
    D = 0.1
    dt = 0.2
    dx = 0.5
    T = 10
    n = 10
    simulation = ClassicFick(dt, dx, T, n, D, lower_bound, upper_bound)
    a = D * dt / (dx ** 2)

    # condicoes de contorno:
    C = [[0] * int(simulation.nodes)]

    for i in np.arange(0, T, dt):
        print("solving t={}s".format(i))
        main_diag = [-(1 + 2 * a)] * (simulation.nodes - 2)
        secondary_diag = [a] * (simulation.nodes - 3)
        l_diag = [0] + secondary_diag
        u_diag = secondary_diag + [0]
        d = create_d(C[int(i / dt)][1:-1], a, simulation.lower_bound)
        solution = TDMA_solver(main_diag, l_diag, u_diag, d)
        solution.insert(0, simulation.lower_bound)
        solution.append(simulation.upper_bound)
        C.append(solution)
        print(solution)

    final_solution = pd.DataFrame(C)
    final_solution.to_excel(r"C:\Users\Julia\Documents\tcc\TF\codes\results\classicFick1013_2.xlsx")


if __name__ == "__main__":
    run()
