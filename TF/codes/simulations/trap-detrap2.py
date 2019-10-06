#!/usr/bin/env python
import numpy as np
import pandas as pd
import math

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


    def boundary_conditions(self, lower, upper):
        self.lower_bound = lower
        self.upper_bound = upper


def run():
    #concentracao superficial constante
    lower_bound = 25
    #concentracao no meio
    upper_bound = 0
    #condicoes de teste
    #coeficiente de difusao
    D = 0.1
    #passo no tempo
    dt = 0.2
    dx = 0.5
    T = 10
    n = 10

    simulation = TrapDetrap(dt, dx, T, n, D, lower_bound, upper_bound)

    # condicoes de contorno:
    N_dif = [[0]*int(simulation.nodes)]
    N_trap = [[0]*int(simulation.nodes)]
    N = [[0]*int(simulation.nodes)]

    for j in np.arange(1, T, dt):
        print("solving t={}s".format(i))
        for i in np.arange(1, n, dx):
            gamma = (N_dif[-1][i]*simulation.trap_concentration - N_dif[-1][i]*N_trap[-1][i] - \
                    simulation.host_atom_concentration*N_trap[-1][i]*simulation.exponential_term)*\
                    simulation.K*simulation.dt

        d = create_d(C[int(i/dt)][1:-1], a, simulation.lower_bound)
        solution = TDMA_solver(main_diag, l_diag, u_diag, d)
        solution.insert(0,simulation.lower_bound)
        solution.append(simulation.upper_bound)
        C.append(solution)
        print(solution)

    final_solution = pd.DataFrame(C)
    final_solution.to_excel(r"C:\Users\Julia\Documents\tcc\TF\codes\results\detrap2.xlsx")
