print("methods")

# Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMA_solver(diag, l_diag, u_diag, d):
    """
    :param diag:
    :param l_diag:
    :param u_diag:
    :param d:
    :return:
    """
    n = len(d)  # number of equations
    u_diag[0] = u_diag[0] / diag[0]
    d[0] = d[0] / diag[0]
    for eq in range(1, n):
        print("TDMA_eq:", eq)
        u_diag[eq] = u_diag[eq] / (diag[eq] - l_diag[eq] * u_diag[eq - 1])
        d[eq] = (d[eq] - l_diag[eq] * d[eq - 1]) / (diag[eq] - l_diag[eq] * u_diag[eq - 1])
    xn = [None] * n
    xn[-1] = d[-1]
    for x in range(n - 2, -1, -1):
        xn[x] = (d[x] - u_diag[x] * xn[x + 1])

    return xn
