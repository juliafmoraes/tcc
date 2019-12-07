#!/usr/bin/env python
import argparse

import numpy as np
import pandas as pd
import math
import openpyxl

from simulations.simulation import Simulation
from methods import TDMA_solver, create_d


class ClassicFickGas(Simulation):

    def __init__(self, dt, dx, T, n, lower, upper):

        super().__init__(dt, dx, T, n)
        self.temperature = 718  # K
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_zero = 8.37 * (10 ** (-8))  # m2/s
        self.beta = 0.0001
        self.c_eq = 2.0358 * (10 ** 28)  # m^-3
        self.D_coef = self.D_zero * math.exp(-self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.boundary_conditions(lower, upper)

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper

    def get_surface_conc(self, t):
        return self.c_eq*(1-math.exp(-self.beta*t))

    def run(self, outdir: str):
        D = self.D_coef
        a = D * dt / (dx ** 2)

        # condicoes de contorno:
        C = [[0] * int(self.nodes)]

        print('total iterations', T // dt)
        main_diag = [-(1 + 2 * a)] * (self.nodes - 2)
        secondary_diag = [a] * (self.nodes - 3)
        l_diag = [0] + secondary_diag
        u_diag = secondary_diag + [0]

        for i in range(1, int(T / dt) + 1):

            d = create_d(C[i % 1000 - 1][1:-1], a, self.get_surface_conc(i*self.dt))
            solution = TDMA_solver(main_diag, l_diag, u_diag, d)
            solution.insert(0, self.get_surface_conc(i*self.dt))
            solution.append(self.upper_bound)
            C.append(solution)

            if i % 1000 == 0:
                print("solving t={}s  --> {}".format(i, i / int(T / dt)))
                C = C[-1:]
            if i*dt in [60.0, 600.0, 1800.0] or i*dt % 3600 == 0:
                final_solution = pd.DataFrame(C)
                final_solution.to_excel(
                    r"{}\classicFickGas_{}.xlsx".format(outdir, i*dt),
                        index=False, header=[x for x in range(0, self.nodes)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Classic Fick - Gas boundary')
    parser.add_argument('outdir', type=str, help='Output dir for results')
    args = parser.parse_args()

    upper_bound = 0
    lower_bound = 0
    dt = 0.01  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 22 * 60 * 60  # s (22horas)
    n = 35 * (10 ** (-6))  # micro -> m
    simulation = ClassicFickGas(dt, dx, T, n, lower_bound, upper_bound)

    simulation.run(args.outdir)
