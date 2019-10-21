#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import openpyxl
import xlrd

from simulations.simulation import Simulation


class TrapDetrap(Simulation):
    def __init__(self, dt, dx, T, n, temperature, trap_concentration, host_atom_concentration, detrapping_energy,
                 diffusion_energy, capture_radius, D_zero):
        super().__init__(dt, dx, T, n)

        # values
        self.temperature = temperature  # K
        self.trap_concentration = trap_concentration  # m-3
        self.host_atom_concentration = host_atom_concentration  # m-3
        self.detrapping_energy = detrapping_energy # J
        self.diffusion_energy = diffusion_energy # J
        self.capture_radius = capture_radius  # m
        self.D_zero = D_zero #m2/s

        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_coef = self.D_zero * math.exp(-1 * self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.k_d = math.exp(-1 * self.detrapping_energy / (self.boltzmann_constant * self.temperature))
        self.K = self.capture_radius * self.D_coef
        self.Fo = self.D_coef * self.dt / (self.dx ** 2)
        self.upper_bound = 0

    def run(self, save_point=None):
        if save_point is None:
            # never save
            save_point = self.iterations+1

        N_dif = [[0] * int(self.nodes)]
        N_trap = [[0] * int(self.nodes)]
        # N = [[0] * int(self.nodes)]

        for j in range(self.iterations):
            if j % 10000 == 0:
                print("solving iteration j={}  t={}".format(j, j*self.dt))

            N_dif_solution = [0] * int(self.nodes)
            N_trap_solution = [0] * int(self.nodes)
            # N_solution = [0] * int(self.nodes)

            for i in range(0, self.nodes - 2):
                gamma = ((N_dif[-1][i] * (self.trap_concentration - N_trap[-1][i])) -
                         (self.k_d * self.host_atom_concentration * N_trap[-1][i])
                        ) * self.K * self.dt
                N_trap_solution[i] = gamma + N_trap[-1][i]
                if i == 0:
                    # N_solution[i] = 2 * self.trap_concentration
                    N_dif_solution[i] = 2 * self.trap_concentration - N_trap_solution[i]
                else:
                    N_dif_solution[i] = N_dif[-1][i] + self.Fo * (N_dif[-1][i + 1] - 2 * N_dif[-1][i] + N_dif[-1][i - 1]) - gamma
                    # N_solution[i] = N_trap_solution[i] + N_dif_solution[i]

            N_dif_solution.append(self.upper_bound)
            N_trap_solution.append(self.upper_bound)
            # N_solution.append(self.upper_bound)

            N_dif.append(N_dif_solution)
            N_trap.append(N_trap_solution)
            # N.append(N_solution)

            #save results to use
            if j*dt in [60, 600, 1800] or j*dt % 3600 == 0:
                writer = pd.ExcelWriter(r"C:\Users\Julia\Documents\tcc\TF\codes\results\new\1\trap_detrap{}.xlsx".format(j*dt))
                df_N_dif = pd.DataFrame(N_dif_solution).T
                df_N_trap = pd.DataFrame(N_trap_solution).T

                df_N_dif.to_excel(writer, 'N_dif', index=False, header=[x for x in range(0, self.nodes+1)])
                df_N_trap.to_excel(writer, 'N_trap', index=False, header=[x for x in range(0, self.nodes+1)])

            #save results for further iterations
            if j > 0 and j % save_point == 0:
                print("save_point: {} for j={}  t={}".format(save_point, j, j*dt))
                writer = pd.ExcelWriter(r"C:\Users\Julia\Documents\tcc\TF\codes\results\new\1\save_points\sp_trap_detrap{}.xlsx".format(j * dt))
                df_N_dif_sp = pd.DataFrame(N_dif[-5:])
                df_N_trap_sp = pd.DataFrame(N_trap[-5:])

                df_N_dif_sp.to_excel(writer, 'N_dif', index=False)
                df_N_trap_sp.to_excel(writer, 'N_trap', index=False)
                writer.save()

                N_dif = N_dif[-5:]
                N_trap = N_trap[-5:]


if __name__ == "__main__":
    temperature = 673  # K
    trap_concentration = 1.31 * (10 ** 28)  # m-3
    host_atom_concentration = 7.29 * (10 ** 28)  # m-3
    detrapping_energy = 4.48609 * (10 ** (-20))  # J
    diffusion_energy = 1.7622 * (10 ** (-19))  # J
    capture_radius = 0.38 * (10 ** (-9))  # m
    D_zero = 8.37 * (10 ** (-8))  # m2/s

    dt = 0.0001  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 2 * 60 * 60  # s (2horas)
    n = 20 * (10 ** (-6))  # micro -> m

    simulation = TrapDetrap(dt, dx, T, n, temperature, trap_concentration, host_atom_concentration, detrapping_energy,
                            diffusion_energy, capture_radius, D_zero)
    simulation.run(200000)
