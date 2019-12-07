#!/usr/bin/env python
import argparse

import numpy as np
import pandas as pd
import math
import openpyxl
import xlrd

from simulations.simulation import Simulation


class TrapDetrap(Simulation):
    def __init__(self, dt, dx, T, n, temperature, trap_concentration, host_atom_concentration, detrapping_energy,
                 diffusion_energy, capture_radius, D_zero, c_eq):
        super().__init__(dt, dx, T, n)

        # input values
        self.temperature = temperature  # K
        self.trap_concentration = trap_concentration  # m-3
        self.host_atom_concentration = host_atom_concentration  # m-3
        self.detrapping_energy = detrapping_energy  # J
        self.diffusion_energy = diffusion_energy  # J
        self.capture_radius = capture_radius  # m
        self.D_zero = D_zero  # m2/s
        self.c_eq = c_eq  # m^-3

        # constants and calculated values
        self.beta = 0.0001
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_coef = self.D_zero * math.exp(-1 * self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.k_d = math.exp(-1 * self.detrapping_energy / (self.boltzmann_constant * self.temperature))
        self.K = self.capture_radius * self.D_coef
        self.Fo = self.D_coef * self.dt / (self.dx ** 2)
        self.upper_bound = 0

    def run(self, outdir: str):
        N_dif = [[0] * int(self.nodes)]
        N_trap = [[0] * int(self.nodes)]

        for j in range(self.iterations + 1):

            N_dif_solution = [0] * int(self.nodes)
            N_trap_solution = [0] * int(self.nodes)

            a = self.k_d * self.host_atom_concentration
            b = self.K * self.dt
            for j in range(int(10801 / self.dt), self.iterations + 1):

                N_dif_solution = [0] * int(self.nodes)
                N_trap_solution = [0] * int(self.nodes)

                for i in range(0, self.nodes - 2):
                    gamma = ((N_dif[-1][i] * (self.trap_concentration - N_trap[-1][i])) -
                             (a * N_trap[-1][i])
                             ) * b
                    N_trap_solution[i] = gamma + N_trap[-1][i]
                    if i == 0:
                        N_dif_solution[i] = self.c_eq * (1 - math.exp(-self.beta * j * self.dt)) - N_trap_solution[i]
                    else:
                        N_dif_solution[i] = N_dif[-1][i] + self.Fo * (
                                N_dif[-1][i + 1] - 2 * N_dif[-1][i] + N_dif[-1][i - 1]) - gamma
                    if N_dif_solution[i] < 10 ** (-6):
                        break

            N_dif_solution.append(self.upper_bound)
            N_trap_solution.append(self.upper_bound)

            N_dif.append(N_dif_solution)
            N_trap.append(N_trap_solution)

            #save results
            if j % 100000 == 0:
                print("solving t={}s  --> {}".format(j * dt, j / int(T / dt)))
                N_dif = N_dif[-5:]
                N_trap = N_trap[-5:]
            if j * dt in [60.0, 600.0, 1800.0] or j * dt % 3600 == 0:
                writer = pd.ExcelWriter(
                    r"{}\trapDetrapGas_{}.xlsx".format(outdir, j * dt))
                df_N_dif = pd.DataFrame(N_dif_solution).T
                df_N_trap = pd.DataFrame(N_trap_solution).T

                df_N_dif.to_excel(writer, 'N_dif', index=False, header=[x for x in range(0, self.nodes + 1)])
                df_N_trap.to_excel(writer, 'N_trap', index=False, header=[x for x in range(0, self.nodes + 1)])

                writer.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trapping-Detrapping - Gas Nitriding boundary')
    parser.add_argument('outdir', type=str, help='Output dir for results')
    args = parser.parse_args()

    temperature = 673  # K
    trap_concentration = 1.31 * (10 ** 28)  # m-3
    host_atom_concentration = 7.29 * (10 ** 28)  # m-3
    detrapping_energy = 4.48609 * (10 ** (-20))  # J
    diffusion_energy = 1.7622 * (10 ** (-19))  # J
    capture_radius = 0.38 * (10 ** (-9))  # m
    D_zero = 8.37 * (10 ** (-8))  # m2/s
    c_eq = 2.0358 * (10 ** 28) # m-3

    dt = 0.0001  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 22 * 60 * 60  # s (22horas)
    n = 20 * (10 ** (-6))  # micro -> m

    simulation = TrapDetrap(dt, dx, T, n, temperature, trap_concentration, host_atom_concentration, detrapping_energy,
                            diffusion_energy, capture_radius, D_zero, c_eq)
    simulation.run(args.outdir)
