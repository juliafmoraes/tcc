#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import openpyxl
import xlrd

from simulations.simulation import Simulation


class TrapDetrap(Simulation):
    def __init__(self, dt, dx, T, n):
        super().__init__(dt, dx, T, n)

        # values
        self.temperature = 673  # K
        self.trap_concentration = 1.31 * (10 ** 28)  # m-3
        self.host_atom_concentration = 7.29 * (10 ** 28)  # m-3
        # self.detrapping_energy = 1.45 #eV
        # self.diffusion_energy = 1.1  # eV
        self.detrapping_energy = 4.48609 * (10 ** (-20))  # J
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        # self.boltzmann_constant = 8.617 * (10 ** (-5))  # eV/K
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.capture_radius = 0.38 * (10 ** (-9))  # m
        # self.D_zero = (10 ** (-7)) #m2/s
        self.D_zero = 8.37 * (10 ** (-8)) #m2/s
        self.D_coef = self.D_zero * math.exp(-self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        # self.k_t = self.D_coef
        # self.k_d = self.D_zero * math.exp(-self.detrapping_energy / (self.boltzmann_constant * self.temperature))
        self.K = 4 * math.pi * self.capture_radius * self.D_coef
        self.Fo = self.D_coef * self.dt / (self.dx ** 2)
        self.upper_boundary = 0


def run():
    # condicoes de teste
    # dt = 0.0001
    # dx = 0.5
    # T = 60
    # n = 10

    # condicoes reais
    # dt = 0.00001  # s
    dt = 0.01  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 2 * 60 * 60  # s (2horas)
    n = 20 * (10 ** (-6))  # micro -> m

    simulation = TrapDetrap(dt, dx, T, n)

    N_dif = [[0] * int(simulation.nodes)]
    N_trap = [[0] * int(simulation.nodes)]
    N = [[0] * int(simulation.nodes)]

    #for slices
    # N = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\3\trap_detrap_8.xlsx', sheet_name='N')
    # N_dif = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\3\trap_detrap_8.xlsx', sheet_name='N_dif')
    # N_trap = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\3\trap_detrap_8.xlsx', sheet_name='N_trap')
    #
    # N = N.values.tolist()
    # N_dif = N_dif.values.tolist()
    # N_trap = N_trap.values.tolist()

    gamma_list = [[0] * int(simulation.nodes)]

    for j in range(simulation.iterations + 1):
        if j == 100001:
            break
        print("solving j={}".format(j))

        N_dif_solution = [0] * int(simulation.nodes)
        N_trap_solution = [0] * int(simulation.nodes)
        N_solution = [0] * int(simulation.nodes)
        gamma_for_iteration = [0] * int(simulation.nodes)

        for i in range(0, simulation.nodes - 1):
            gamma = ((N_dif[-1][i] * (simulation.trap_concentration - N_trap[-1][i])) -
                     (N_trap[-1][i] * simulation.D_coef * simulation.host_atom_concentration)
                     ) * simulation.K * simulation.dt/1000
            gamma_for_iteration[i] = gamma
            N_trap_solution[i] = gamma + N_trap[-1][i]

            if i == 0:
                N_solution[i] = 2 * 1.31 * (10 ** 28)  #concentracao na superficie cte
                N_dif_solution[i] = N_solution[i] - N_trap_solution[i]
            else:
                N_dif_solution[i] = N_dif[-1][i] + simulation.Fo * (N_dif[-1][i + 1] - 2 * N_dif[-1][i] + N_dif[-1][i - 1]) - gamma
                N_solution[i] = N_trap_solution[i] + N_dif_solution[i]

        N_dif_solution.append(simulation.upper_boundary)
        N_trap_solution.append(simulation.upper_boundary)
        N_solution.append(simulation.upper_boundary)

        N_dif.append(N_dif_solution)
        N_trap.append(N_trap_solution)
        N.append(N_solution)
        gamma_list.append(gamma_for_iteration)

    # final_N_solution = pd.DataFrame(N)
    # final_N_dif_solution = pd.DataFrame(N_dif)
    # final_N_trap_solution = pd.DataFrame(N_trap)

    # final_solution_excel = final_solution.iloc[[x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]
    # final_solution_excel.to_excel(r"C:\Users\Julia\Documents\tcc\TF\codes\results\trap_detrap2_2.xlsx")

    # final_N_solution_excel = final_N_solution.iloc[
    #     [x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]
    # final_N_dif_solution_excel = final_N_dif_solution.iloc[
    #     [x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]
    # final_N_trap_solution_excel = final_N_trap_solution.iloc[
    #     [x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]

    # final_N_solution_excel = final_N_solution.iloc[
    #     [x for x in range(0, len(N), int(len(N) / 50))]]
    # final_N_dif_solution_excel = final_N_dif_solution.iloc[
    #     [x for x in range(0, len(N), int(len(N) / 50))]]
    # final_N_trap_solution_excel = final_N_trap_solution.iloc[
    #     [x for x in range(0, len(N), int(len(N) / 50))]]

#for slices
    final_N_solution = pd.DataFrame(N[-5:])
    final_N_dif_solution = pd.DataFrame(N_dif[-5:])
    final_N_trap_solution = pd.DataFrame(N_trap[-5:])

    final_N_solution_excel = final_N_solution
    final_N_dif_solution_excel = final_N_dif_solution
    final_N_trap_solution_excel = final_N_trap_solution

    writer = pd.ExcelWriter(r"C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\4\trap_detrap_1.xlsx")
    final_N_solution_excel.to_excel(writer, 'N', index=False)
    final_N_dif_solution_excel.to_excel(writer, 'N_dif', index=False)
    final_N_trap_solution_excel.to_excel(writer, 'N_trap', index=False)
    writer.save()


if __name__ == "__main__":
    run()
