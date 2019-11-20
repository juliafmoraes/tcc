# Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMA_solver(diag, l_diag, u_diag, d):
    """
    :param diag: main diagonal
    :param l_diag: lower diag
    :param u_diag: upper
    :param d: right side
    :return:
    """
    u_diag_copy = u_diag.copy()
    d_copy = d.copy()
    n = len(d)  # number of equations
    u_diag_copy[0] = u_diag_copy[0] / diag[0]
    d_copy[0] = d[0] / diag[0]
    for eq in range(1, n - 1):
        u_diag_copy[eq] = u_diag_copy[eq] / (diag[eq] - l_diag[eq] * u_diag_copy[eq - 1])
        d_copy[eq] = (d_copy[eq] - l_diag[eq] * d_copy[eq - 1]) / (diag[eq] - l_diag[eq] * u_diag_copy[eq - 1])

    eq = n - 1
    d_copy[eq] = (d_copy[eq] - l_diag[eq] * d_copy[eq - 1]) / (diag[eq] - l_diag[eq] * u_diag_copy[eq - 1])
    xn = [None] * n
    xn[-1] = d_copy[-1]
    for x in range(n - 2, -1, -1):
        xn[x] = (d_copy[x] - u_diag_copy[x] * xn[x + 1])

    return xn


def create_d(c, a, lb):
    d = [-1 * element for element in c]
    d[0] = d[0] - a * lb
    return d
