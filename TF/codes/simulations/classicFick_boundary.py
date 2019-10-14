#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import openpyxl

from simulations.simulation import Simulation
from methods import TDMA_solver


class ClassicFickBoundary(Simulation):

    def __init__(self, dt, dx, T, n, lower, upper):
        """
        :type D: float
        """
        super().__init__(dt, dx, T, n)
        self.temperature = 673  # K
        self.alfa = 1
        self.j_zero = 3 * (10**(-9)) #m/2?
        # self.host_atom_concentration = 50
        self.host_atom_concentration = 7.29 * (10 ** 28)  # m-3
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_zero = 8.37 * (10 ** (-8))  # m2/s
        self.D_coef = self.D_zero * math.exp(-self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        # self.j_zero = 1 / 2
        self.flux_term = self.alfa * self.j_zero * self.dt
        self.boundary_conditions(lower, upper)

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper


    def create_d(self, c, a, lb):
        d = [-1 * element for element in c]
        d[0] = c[0] + self.flux_term*self.host_atom_concentration
        d[1] = -a*(c[0] + self.flux_term * self.host_atom_concentration) - c[1]
        return d


def run():
    # condicoes de teste
    lower_bound = 0
    upper_bound = 0

    # D = 0.1
    # dt = 0.2
    # dx = 0.5
    # T = 10
    # n = 10
    # condicoes reais
    # dt = 0.00001  # s
    dt = 0.01  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 2 * 60 * 60  # s (2horas)
    n = 20 * (10 ** (-6))  # micro -> m
    simulation = ClassicFickBoundary(dt, dx, T, n, lower_bound, upper_bound)
    D = simulation.D_coef
    a = D * dt / (dx ** 2)

    # condicoes de contorno:
    C = [[0] * int(simulation.nodes)]

    for i in np.arange(0, T, dt):

        print("solving t={}s".format(i))
        main_diag = [1 - a + simulation.flux_term, -1*(a**2+2*a+1)] + [-(1 + 2 * a)] * (simulation.nodes - 3)
        secondary_diag = [a] * (simulation.nodes - 3)

        u_diag = [a] + secondary_diag + [0]
        l_diag = [0, a ** 2 - a*simulation.flux_term] + secondary_diag

        d = simulation.create_d(C[int(i / dt)][0:-1], a, simulation.lower_bound)
        solution = TDMA_solver(main_diag, l_diag, u_diag, d)
        # solution.insert(0, simulation.lower_bound)
        solution.append(simulation.upper_bound)
        C.append(solution)
        print(solution)

    final_solution = pd.DataFrame(C)
    final_solution.to_excel(r"C:\Users\Julia\Documents\tcc\TF\codes\results\classicFick_not_constant_real.xlsx")


if __name__ == "__main__":
    run()
