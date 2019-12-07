#!/usr/bin/env python

import pandas as pd


class Simulation:

    def __init__(self, dt, dx, T, n):
        self.dt = dt
        self.dx = dx
        self.T = T
        self.n = n
        self.nodes: int = int(n / dx)
        self.iterations: int = int(T / dt)
        self.results = pd.DataFrame()


if __name__ == "__main__":
    pass
