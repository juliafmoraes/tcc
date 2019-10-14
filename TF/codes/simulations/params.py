#!/usr/bin/env python

import pandas as pd
from simulations.simulation import Simulation
import math

class Parameters(Simulation):

    def __init__(self, dt, dx, T, n):
        super().__init__(dt, dx, T, n)

        self.temperature = 673  # K
        self.host_atom_concentration = 7.29 * (10 ** 28)  # m-3
        self.diffusion_energy = 1.7622 * (10 ** (-19))  # J
        self.boltzmann_constant = 1.38064852 * (10 ** (-23))  # m2 kg s-2 K-1 = (J/K)
        self.D_zero = 8.37 * (10 ** (-8))  # m2/s
        self.D_coef = self.D_zero * math.exp(-self.diffusion_energy / (self.boltzmann_constant * self.temperature))
        self.Fo = self.D_coef * self.dt / (self.dx ** 2)
        self.upper_boundary = 0

if __name__ == "__main__":
    pass
