#!/usr/bin/env python
import numpy as np
import pandas as pd
import math

from simulations.simulation import Simulation
from methods import TDMA_solver


class TrapDetrap(Simulation):
    def __init__(self, dt, dx, T, n, D: float, lower, upper):
        """
        :type D: float
        """
        super().__init__(dt, dx, T, n)
        self.D_coef = D
        self.boundary_conditions(lower, upper)

        self.available = [[0] * int(self.nodes)]
        self.trapped = [[0] * int(self.nodes)]
        self.total = self.get_total_concentration()

        self.cr_concentration = 19.5
        self.exp_detrap = math.exp((-0.28 * 1.60218 * 10 ** -19) / ((1.38 * 10 ** -23) * 673))

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper

    def get_total_concentration(self):
        return [map(sum, zip(*t)) for t in zip(self.available, self.trapped)]

    def diffuse(self, i):
        C = self.available
        a = self.D_coef * self.dt / (self.dx ** 2)

        print("solving t={}s".format(i))
        main_diag = [-(1 + 2 * a)] * (self.nodes - 2)
        secondary_diag = [a] * (self.nodes - 3)
        u_diag = [0]
        u_diag.extend(secondary_diag)
        l_diag = secondary_diag
        l_diag.append(0)
        d = create_d(C[int(i / self.dt)][1:-1], a, self.lower_bound)
        solution = TDMA_solver(main_diag, l_diag, u_diag, d)
        solution.insert(0, self.lower_bound)
        solution.append(self.upper_bound)

        C.append(solution)
        print(solution)

        self.available = C

    def trap(self):
        T = self.trapped[-1]
        for i, node in enumerate(self.available[-1]):
            if node > self.cr_concentration:
                self.available[-1][i] = node - self.cr_concentration
                T[i] = self.cr_concentration
        self.trapped.append(T)

    def detrap(self):
        for i, node in enumerate(self.trapped[-1]):
            if node > 0:
                self.available[-1][i] += node * self.exp_detrap
                self.trapped[-1][i] -= node * self.exp_detrap


def create_d(c, a, lb):
    d = [-1 * element for element in c]
    d[0] = d[0] - a * lb
    return d


def run():
    lower_bound = 25
    upper_bound = 0
    D = 0.1
    dt = 0.2
    dx = 0.05
    T = 3600*10
    n = 10
    simulation = TrapDetrap(dt, dx, T, n, D, lower_bound, upper_bound)
    a = D * dt / (dx ** 2)

    # condicoes de contorno:
    C = [[0] * int(simulation.nodes)]
    simulation.available = C
    for i in np.arange(0, T, dt):
        simulation.diffuse(i)
        simulation.trap()
        # simulation.detrap()

    final_solution = pd.DataFrame(simulation.get_total_concentration())
    final_solution_excel = final_solution.iloc[[x for x in range(0,simulation.iterations, int(simulation.iterations/50))]]
    final_solution_excel.to_excel(r"C:\Users\Julia\Documents\tcc\TF\codes\results\trap-detrap.xlsx")


if __name__ == "__main__":
    run()
