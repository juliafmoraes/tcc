#!/usr/bin/env python

def from_vol_to_at(input_N):
    NA =  6.0221409 * (10**23)

    #composition %
    comp_Ni = 0.113
    comp_Cr = 0.195
    comp_Fe = 1 - comp_Ni - comp_Cr

    #molar volumes m3/mol
    V_Ni = 6.59 * (10 ** (-6))
    V_Cr = 7.23 * (10 ** (-6))
    V_Fe = 7.09 * (10 ** (-6))

    TOTAL = NA/(comp_Ni*V_Ni + comp_Cr*V_Cr + comp_Fe*V_Fe) #at/m3s

    return input_N / (input_N + TOTAL)

# print(from_vol_to_at(2.62*(10**28)))