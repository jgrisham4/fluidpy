#!/usr/bin/env python

import sys
import os
import chemistry

################################################################
# Inputs
################################################################

# Defining temperature and pressure
p_m = 13.0    # atm
T = 4500.0    # K
#T = 3000.0    # K

MW_N = 14.0067  # g/mol or apparently kg/(kg-mol)
MW_O = 15.9994  # kg/(kg-mol)
deltaHf_N = 4.714e8
deltaHf_O = 2.47e8
deltahf_N = deltaHf_N/MW_N
deltahf_O = deltaHf_O/MW_O

#                      _________________________________
# Defining species    |____MW____|___h_f____|_name_|_nu_|
N2 = chemistry.species(28.0134e-3,       0.0,  "N2", 1.0, "diatomic", theta_v=3390.0)
O2 = chemistry.species(31.9988e-3,       0.0,  "O2", 1.0, "diatomic", theta_v=2270.0)
N  = chemistry.species(14.0067e-3, deltahf_N,   "N", 2.0, "monatomic")
O  = chemistry.species(15.9994e-3, deltahf_O,   "O", 2.0, "monatomic")

# Defining reactions
O2_dissociation = chemistry.reaction([O2], [O], "O2 -> 2O",C_c=1.2e6, 
    eta_c=-0.5, Theta=59500.0)
N2_dissociation = chemistry.reaction([N2], [N], "N2 -> 2N",C_c=1.8e4, 
    eta_c=0.0 , Theta=113000.0)

# Defining reacting mixture
mixture = chemistry.air([O2_dissociation, N2_dissociation], tol=1.0e-10)

# Finding equilibrium at given pressure and temperature
#mixture.find_equilibrium(0.1*101.325e3, 4500.0)
mixture.find_equilibrium(p_m, T)

# Writing history to file
with open("hw3prob3.dat", "w") as outfile:
    outfile.write("{0:5} {1:10} {2:10} {3:10} {4:10}\n"
        .format("Iter", "p_O", "p_N", "p_O2", "p_N2"))
    for n in range(0, len(mixture.hist[0])):
        outfile.write("{0:5} {1:1.4e} {2:1.4e} {3:1.4e} {4:1.4e}\n"
            .format(n, mixture.hist[0][n], mixture.hist[1][n], 
              mixture.hist[2][n], mixture.hist[3][n]))

# Performing a check
p_total = 0.0
for ph in mixture.hist:
    p_total += ph[-1]
print("\nCheck using Dalton's law: 0.1 ?= {:6f}\n".format(p_total))

# Writing results to file
with open("hw3prob3_pp.tex", "w") as fp:
    fp.write("\\begin{align*}\n")
    fp.write("p_\\text{O} &= ")
    fp.write("{:1.4f}".format(O.p))
    fp.write(" \\ \\text{atm}\\\\\n")
    fp.write("p_\\text{N} &= ")
    fp.write("{:1.4f}".format(N.p))
    fp.write(" \\ \\text{atm}\\\\\n")
    fp.write("p_\\text{O$_2$} &= ")
    fp.write("{:1.6f}".format(O2.p))
    fp.write(" \\ \\text{atm}\\\\\n")
    fp.write("p_\\text{N$_2$} &= ")
    fp.write("{:1.4f}".format(N2.p))
    fp.write(" \\ \\text{atm}\n")
    fp.write("\\end{align*}\n")

# Finding properties for each species
species_list = [O, N, O2, N2]
MW_mixture = 0.0
h_m = 0.0
e_m = 0.0
for s in species_list:
    s.X = s.p/p_m
    MW_mixture += s.X*s.MW
for s in species_list:
    s.c = s.X*s.MW/MW_mixture
    current_h = s.h(T)
    current_e = s.e(T)
    h_m += s.c*current_h
    e_m += s.c*current_e

print("h_m({}) = {:1.3e} J/kg".format(T, h_m))
print("e_m({}) = {:1.3e} J/kg".format(T, e_m))

# Writing results to file
#with open("hw3prob3_eh.tex", "w") as fp:
#    fp.write("\\begin{equation}\\boxed{\n")
#    fp.write("e_m(")
#    fp.write("{:4.0f}) = \sum c_i e_i = ".format(T))
#    fp.write("{:1.3f} \\ \\text".format(e_m/1.0e6))
#    fp.write("{MJ/kg}\n")
#    fp.write("}\\end{equation}\n")
#    fp.write("\\begin{equation}\\boxed{\n")
#    fp.write("h_m(")
#    fp.write("{:4.0f}) = \sum c_i h_i = ".format(T))
#    fp.write("{:1.3f} \\ \\text".format(h_m/1.0e6))
#    fp.write("{MJ/kg}\n")
#    fp.write("}\\end{equation}\n")
#    fp.write("where $h_i = e_i + R_i T$.\n")
