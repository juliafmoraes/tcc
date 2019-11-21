#!/usr/bin/env python
import argparse

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
        self.j = 3 * (10 ** (-9))  # m/2?
        self.q = 1.602 * (10 ** (-19))
        self.host_surface_concentration = 8 * (10 ** 25)  # at/m2
        # self.host_atom_concentration = 50
        self.host_atom_concentration = 7.29 * (10 ** 28)  # m-3
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_zero = 8.37 * (10 ** (-8))  # m2/s
        self.D_coef = self.D_zero * math.exp(-self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.j_zero = self.j / (self.q * self.host_surface_concentration)
        self.flux_term = self.alfa * self.j_zero * self.dt
        self.boundary_conditions(lower, upper)

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper

    def create_d(self, c):
        d = [element for element in c]
        d[0] = d[0] + self.flux_term * self.host_atom_concentration
        return d

    def run(self, outdir: str):

        D = self.D_coef
        a = D * dt / (dx ** 2)

        # condicoes de contorno:
        C = [[0] * int(self.nodes)]

        print('total iterations', T // dt)
        main_diag = [1 + a + self.flux_term] + [(1 + 2 * a)] * (self.nodes - 2)
        secondary_diag = [-a] * (self.nodes - 2)

        u_diag = secondary_diag + [0]
        l_diag = [0] + secondary_diag
        for i in range(1, self.iterations + 1):
            d = self.create_d(C[i % 1000 - 1][0:-1])  # %1000 because we only maintain a 1000 len sized list
            solution = TDMA_solver(main_diag, l_diag, u_diag, d)
            solution.append(self.upper_bound)
            C.append(solution)

            if i % 1000 == 0:
                print("solving t={}s  --> {}".format(i, i / int(T / dt)))
                C = C[-1:]
            if i % 10000 == 0:
                final_solution = pd.DataFrame(C)
                final_solution.to_excel(
                    r"{}\classicFickPlasma_teste{}.xlsx".format(outdir, i),
                    index=False, header=[x for x in range(0, self.nodes)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Classic Fick - Plasma boundary')
    parser.add_argument('outdir', type=str, help='Output dir for results')
    args = parser.parse_args()

    lower_bound = None
    upper_bound = 0

    dt = 0.01  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 22 * 60 * 60  # s (22horas)
    n = 20 * (10 ** (-6))  # micro -> m
    simulation = ClassicFickBoundary(dt, dx, T, n, lower_bound, upper_bound)

    simulation.run(args.outdir)
