#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import openpyxl

from simulations.simulation import Simulation


class TrapDetrap(Simulation):
    def __init__(self, dt, dx, T, n, D: float, lower, upper):
        super().__init__(dt, dx, T, n)
        self.D_coef = D
        self.boundary_conditions(lower, upper)

        # values
        self.temperature = 400 + 273    #K
        self.trap_concentration = 1.31 * (10**22)   #cm-3
        self.host_atom_concentration = 7.29 * (10**22)  #cm-3
        self.detrapping_energy = 4.48609 * (10**(-18))   #J
        self.diffusion_energy = 1.7622 * (10**(-19)) #J
        self.boltzmann_constant = 1.38064852 * (10**(-23)) #m2 kg s-2 K-1
        self.capture_radius = 0.38 * (10**(-7))
        self.K = 4*math.pi*self.capture_radius*self.D_coef
        self.exponential_term = math.exp(-self.detrapping_energy/(self.boltzmann_constant*self.temperature))
        self.Fo = self.D_coef*self.dt/(self.dx**2)

    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper


def run():
    #concentracao superficial constante
    lower_bound = 3*7.29*(10**(-7+22))
    #concentracao no meio
    upper_bound = 0
    #condicoes de teste
    #coeficiente de difusao - teste
    D = 0.1
    #coeficiente de difusao - real
    # D = 4.81 * (10**(-12))
    #passo no tempo
    dt = 0.2
    dx = 0.5
    T = 10
    n = 10

    simulation = TrapDetrap(dt, dx, T, n, D, lower_bound, upper_bound)

    # condicoes de contorno:
    N_dif =[[simulation.lower_bound]+[0]*int(simulation.nodes-1)]
    N_trap = [[0]*int(simulation.nodes)]
    N = [[simulation.lower_bound]+[0]*int(simulation.nodes-1)]

    for j in range(simulation.iterations):
        print("solving t={}s".format(j))
        N_dif_solution = [simulation.lower_bound]+[0]*int(simulation.nodes-1)
        N_trap_solution = [0]*int(simulation.nodes)
        N_solution = [simulation.lower_bound]+[0]*int(simulation.nodes-1)
        for i in range(1, simulation.nodes-1):
            gamma = (N_dif[-1][i] - N_dif[-1][i]*N_trap[-1][i]/simulation.trap_concentration -
                     simulation.host_atom_concentration*N_trap[-1][i]*simulation.exponential_term)*\
                    simulation.K*simulation.dt
            N_trap_solution[i] = gamma + N_trap[-1][i]
            N_dif_solution[i] = simulation.Fo*(N_dif[-1][i+1]-2*N_dif[-1][i]+N_dif[-1][i-1]) +\
                                N_dif[-1][i] - simulation.dt*gamma
            N_solution[i] = N[-1][i] - N_dif[-1][i] - N_trap[-1][i] + N_trap_solution[i] + N_dif_solution[i]

        N_dif_solution.append(simulation.upper_bound)
        N_trap_solution.append(simulation.upper_bound)
        N_solution.append(simulation.upper_bound)

        N_dif.append(N_dif_solution)
        N_trap.append(N_trap_solution)
        N.append(N_solution)

        # print(solution)

    final_solution = pd.DataFrame(N)
    final_solution.to_excel(r"C:\Users\Ju\tcc\TF\codes\results\trap_detrap2.xlsx")


if __name__ == "__main__":
    run()