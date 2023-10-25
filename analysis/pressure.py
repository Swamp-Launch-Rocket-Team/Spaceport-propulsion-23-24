import numpy as np
import math
import matplotlib.pyplot as plt

from analysis import parameters as par

def chamber_pressure():
    # this function should call the startup(), steady_state(), and tail_off() functions
    None

def startup():

    # Use ideal gas to find gas density?

    None


def steady_state(params: par.Params(), dt: float, dx: float):
    Ab = find_Ab(params, dx)
    C = params.prop.burn_coef * params.prop.dens / math.sqrt(params.prop.heat_ratio / params.prop.R / params.prop.burn_temp * (2 / (params.prop.heat_ratio + 1))**((params.prop.heat_ratio + 1)/(params.prop.heat_ratio - 1)))
    P = (Ab / params.nozzle.throat_area * C)**(1/(1 - params.prop.burn_exp))

    return P


def tail_off():



    None

# Add other functions down here

def find_Ab(params: par.Params(), dt: float) -> np.array:
    
    # diam_burn = (params.prop.Do - params.prop.Di) / 2 / dx
    # length_burn = params.prop.L / 2 / dx
    # N = round(diam_burn) if diam_burn < length_burn else round(length_burn)
    # Ab = np.zeros(N+1, dtype=float)
    # for i in range(N+1):
    #     delta_L = params.prop.burn_rate * params.dt * i
    #     Di = params.prop.Di + 2 * dx * i
    #     L = params.prop.L - 2 * dx * i
    #     Ab[i] = np.pi / 2 * (params.prop.Do**2 - Di**2) + np.pi * Di * L

    diam_burn = (params.prop.Do - params.prop.Di) / 2 / dt
    length_burn = params.prop.L / 2 / dt

    t_diam = (params.prop.Do - params.prop.Di) / 2 / params.prop.burn_rate
    t_length = params.prop.L / 2 / params.prop.burn_rate

    params.dt = dt

    t_burn = t_diam if t_diam < t_length else t_length

    N = round(t_burn / dt)

    print(f"t: {t_burn}, N: {N}")

    Ab = np.zeros(N+1, dtype=float)
    for i in range(N+1):
        delta_L = params.prop.burn_rate * params.dt * i
        Di = params.prop.Di + 2 * delta_L * i
        L = params.prop.L - 2 * delta_L * i
        Ab[i] = np.pi / 2 * (params.prop.Do**2 - Di**2) + np.pi * Di * L

    return Ab

