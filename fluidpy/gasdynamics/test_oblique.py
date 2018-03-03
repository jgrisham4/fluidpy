"""
This script plots the theta-beta-Mach diagram for an oblique shock.
"""

import obliqueshock
import numpy as np
import matplotlib.pyplot as plt

machach = [1.4, 1.6, 2.0, 3.0]
theta = np.linspace(0.0, 60.0, 500)
beta_vec = [[], [], [], []]
theta_vec = [[], [], [], []]
ctr = 0
for mach in machach:
    for th in theta:
        btmp = obliqueshock.wave_angle(mach, th*np.pi/180.0)
        if btmp.imag == 0.0:
            beta_vec[ctr].append(btmp*180.0/np.pi)
            theta_vec[ctr].append(th)
    ctr += 1

fig = plt.figure()
for i in range(0, 3):
    plt.plot(beta_vec[i], theta_vec[i], label="mach={}".format(machach[i]))
plt.legend(loc=1)
plt.grid()


mach = 3.0
thetaprime = 12.0
p2qp1 = obliqueshock.pressure_ratio(mach, thetaprime*np.pi/180.0)
rho2qrho1 = obliqueshock.density_ratio(mach, thetaprime*np.pi/180.0)
t2qt1 = obliqueshock.temperature_ratio(mach, thetaprime*np.pi/180.0)

print("For mach = {} and theta = {}".format(mach, thetaprime))
print("p2/p1 = {0:2.4f} rho2/rho1 = {1:2.4f} T2/T1 = {2:2.4f}".format(p2qp1, rho2qrho1, t2qt1))

plt.show()
