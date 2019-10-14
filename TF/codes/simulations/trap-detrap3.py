#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import openpyxl

from simulations.simulation import Simulation


class TrapDetrap(Simulation):
    def __init__(self, dt, dx, T, n):
        super().__init__(dt, dx, T, n)

        # values
        self.temperature = 673  # K
        self.trap_concentration = 1.31 * (10 ** 28)  # m-3
        self.host_atom_concentration = 7.29 * (10 ** 28)  # m-3
        self.detrapping_energy = 4.48609 * (10 ** (-20))  # J
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.capture_radius = 0.38 * (10 ** (-9))  # m
        self.D_zero = 8.37 * (10 ** (-8))
        self.D_coef = self.D_zero * math.exp(-1 * self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.K = 4 * math.pi * self.capture_radius * self.D_coef
        self.exponential_term = math.exp(-1 * self.detrapping_energy / (self.boltzmann_constant * self.temperature))
        self.Fo = self.D_coef * self.dt / (self.dx ** 2)
        self.j_zero = 3 * (10 ** (-7))  # ??
        self.upper_bound = 0


def run():
    # condicoes de teste
    # coeficiente de difusao
    # D = 0.1
    # #passo no tempo
    # dt = 0.2
    # dx = 0.5
    # T = 10
    # n = 10

    # condicoes reais
    # passo no tempo
    dt = 0.001  # s
    dx = 0.1 * (10 ** (-6))  # micro -> m
    T = 2 * 60 * 60  # s (2horas)
    n = 20 * (10 ** (-6))  # micro -> m

    # dt = 0.001  # s
    # dx = 0.1 * (10 ** (-6))  # micro -> m
    # T = 2 * 60 * 60  # s (2horas)
    # n = 20 * (10 ** (-6))  # micro -> m

    simulation = TrapDetrap(dt, dx, T, n)

    N_dif = [[0] * int(simulation.nodes)]
    N_trap = [[0] * int(simulation.nodes)]
    N = [[0] * int(simulation.nodes)]

    for j in range(simulation.iterations):
        print("solving j={}".format(j))

        N_dif_solution = [0] * int(simulation.nodes)
        N_trap_solution = [0] * int(simulation.nodes)
        N_solution = [0] * int(simulation.nodes)

        for i in range(0, simulation.nodes - 1):
            gamma = (
                            (N_dif[-1][i] * (simulation.trap_concentration - N_trap[-1][i])) -
                            (simulation.host_atom_concentration * N_trap[-1][i] * simulation.exponential_term)
                    ) * simulation.K * simulation.dt
            N_trap_solution[i] = gamma + N_trap[-1][i]
            if i == 0:
                N_dif_solution[i] = -simulation.Fo * (N_dif[-1][i + 1] - N_dif[-1][i]) - gamma + \
                                    simulation.j_zero * simulation.dt * (
                                                simulation.host_atom_concentration - N_dif[-1][i] - N_trap[-1][i])
            else:
                N_dif_solution[i] = simulation.Fo * (N_dif[-1][i + 1] - 2 * N_dif[-1][i] + N_dif[-1][i - 1]) + \
                                    N_dif[-1][i] - simulation.dt * gamma

        N_solution[i] = N_trap_solution[i] + N_dif_solution[i]

        N_dif_solution.append(simulation.upper_bound)
        N_trap_solution.append(simulation.upper_bound)
        N_solution.append(simulation.upper_bound)

        N_dif.append(N_dif_solution)
        N_trap.append(N_trap_solution)
        N.append(N_solution)

    final_N_solution = pd.DataFrame(N)
    final_N_dif_solution = pd.DataFrame(N_dif)
    final_N_trap_solution = pd.DataFrame(N_trap)

    # final_solution_excel = final_solution.iloc[[x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]
    # final_solution_excel.to_excel(r"C:\Users\Julia\Documents\tcc\TF\codes\results\trap_detrap2_2.xlsx")

    final_N_solution_excel = final_N_solution.iloc[
        [x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]
    final_N_dif_solution_excel = final_N_dif_solution.iloc[
        [x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]
    final_N_trap_solution_excel = final_N_trap_solution.iloc[
        [x for x in range(0, simulation.iterations, int(simulation.iterations / 50))]]

    writer = pd.ExcelWriter(r"C:\Users\Julia\Documents\tcc\TF\codes\results\trap_detrap3.xlsx")
    final_N_solution_excel.to_excel(writer, 'N')
    final_N_dif_solution_excel.to_excel(writer, 'N_dif')
    final_N_trap_solution_excel.to_excel(writer, 'N_trap')
    writer.save()


if __name__ == "__main__":
    run()
