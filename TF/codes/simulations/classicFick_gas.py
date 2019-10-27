#!/usr/bin/env python
import argparse

import numpy as np
import pandas as pd
import math
import openpyxl

from simulations.simulation import Simulation
from methods import TDMA_solver


class ClassicFickGas(Simulation):

    def __init__(self, dt, dx, T, n, lower, upper):

        super().__init__(dt, dx, T, n)
        self.temperature = 718  # K
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_zero = 8.37 * (10 ** (-8))  # m2/s
        self.beta = 0.0001
        self.c_eq = 2.8496 * (10**28)
        self.D_coef = self.D_zero * math.exp(-self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.boundary_conditions(lower, upper)

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper

    def get_surface_conc(self, t):
        return self.c_eq*(1-math.exp(-self.beta*t))


def create_d(c, a, simulation, t):
    d = [-1 * element for element in c]
    d[0] = d[0] - a * simulation.get_surface_conc(t)
    return d


def run(outdir: str):
    # condicoes de teste
    lower_bound = 0
    upper_bound = 0
    # D = 0.1
    # dt = 0.2
    # dx = 0.5
    # T = 10
    # n = 10
    # condicoes reais
    dt = 0.01  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 2 * 60 * 60  # s (2horas)
    n = 20 * (10 ** (-6))  # micro -> m
    simulation = ClassicFickGas(dt, dx, T, n, lower_bound, upper_bound)
    D = simulation.D_coef
    a = D * dt / (dx ** 2)

    # condicoes de contorno:
    C = [[0] * int(simulation.nodes)]

    print('total iterations', T // dt)
    main_diag = [-(1 + 2 * a)] * (simulation.nodes - 2)
    secondary_diag = [a] * (simulation.nodes - 3)
    l_diag = [0] + secondary_diag
    u_diag = secondary_diag + [0]

    for i in range(1, int(T / dt) + 1):
        print("solving t={}s  --> {}".format(i, i/int(T / dt)))
        d = create_d(C[i % 1000 - 1][1:-1], a, simulation, i*simulation.dt)
        solution = TDMA_solver(main_diag, l_diag, u_diag, d)
        solution.insert(0, simulation.get_surface_conc(i*simulation.dt))
        solution.append(simulation.upper_bound)
        C.append(solution)
        # print(solution)

        if i % 1000 == 0:
            C = C[-1:]
        if i*dt in [60.0, 600.0, 1800.0] or i*dt % 3600 == 0:
            final_solution = pd.DataFrame(C)
            final_solution.to_excel(
                r"{}\classicFickGas_{}.xlsx".format(outdir, i*dt),
                    index=False, header=[x for x in range(0, simulation.nodes)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Classic Fick - Gas boundary')
    parser.add_argument('outdir', type=str, help='Output dir for results')
    args = parser.parse_args()
    run(args.outdir)
