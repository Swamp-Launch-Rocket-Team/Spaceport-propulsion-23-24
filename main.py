import numpy as np
import matplotlib.pyplot as plt

from analysis import parameters as par
from analysis import nozzle as na
from analysis import pressure as pa
from analysis import ht as ht
from analysis import combustion as cb

# Main code for analysis will go here

#params = par.Params()
# params.prop.heat_ratio = 1.1317
# params.prop.burn_temp = 1583.066
# params.prop.dens = 1.811 * 1000
# params.prop.burn_coef = 0.000226
# params.prop.burn_exp = 0.424
# params.prop.KN = 180
# params.prop.L = 0.5
# params.prop.Do = 0.0508
# params.prop.Di = 0.0254
# params.prop.R = 226.979
# params.prop.burn_rate = params.prop.burn_coef * 3447000 ** params.prop.burn_exp
# params.prop.burn_rate = 0.005
# params.nozzle.throat_area = np.pi * 0.01**2 / 4

# print(params.prop.burn_rate)

# Ab = pa.find_Ab(params, 0.001)
# # print(Ab)
# P = pa.steady_state(params, 100, 0.001)
# plt.figure(1)
# plt.plot(Ab)
# plt.grid()
# plt.show(block=False)
# # plt.pause(0.001)

# plt.figure(2)
# plt.plot(P)
# plt.grid()
# plt.show()
# None

[t,p0] = pa.chamber_pressure()

plt.figure(1)
plt.plot(t[0:-1],p0)
plt.grid()
plt.show()
None