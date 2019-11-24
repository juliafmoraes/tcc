#!/usr/bin/env python
NA = 6.0221409 * (10 ** 23)


def from_vol_to_at(input_N):
    # composition %
    comp_Ni = 0.113
    comp_Cr = 0.195
    comp_Fe = 1 - comp_Ni - comp_Cr

    # molar volumes m3/mol
    V_Ni = 6.59 * (10 ** (-6))
    V_Cr = 7.23 * (10 ** (-6))
    V_Fe = 7.09 * (10 ** (-6))

    TOTAL = NA / (comp_Ni * V_Ni + comp_Cr * V_Cr + comp_Fe * V_Fe)  # at/m3s

    return input_N / (input_N + TOTAL)


# print(from_vol_to_at(2.62*(10**28)))
# print(from_vol_to_at(3.31*(10**28)))

def from_yn_to_m3(yn):
    return (4 * yn) / (2.8147 * (10 ** (-29)) * yn + 4.7134 * (10 ** (-29)))


# print(from_yn_to_m3(0.28))
