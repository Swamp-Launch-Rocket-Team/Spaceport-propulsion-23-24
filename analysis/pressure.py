import numpy as np
import math
import matplotlib.pyplot as plt

from analysis import parameters as par



#returns an array of absolute chamber pressure (MPa) and time (seconds) for use in constructing a graph
# def chamber_pressure(params: par.Params() )-> np.array:
def chamber_pressure():
    # euler's method to find pressure based off of Nakka srm excel spreadsheet for rocket motor design

    #variables
    patm = 0.101 #local atmospheric pressure (MPa)
    erate = 0 #erosion rate (mm/s)
    dc = 75 #chamber diameter (inside) (mm)   
    lc = 470 #chamber length (inside) (mm)
    gstar = 6 #Ratio of core to throat cross-sectional areas, above which no erosive burning occurs (unitless)
    a = 17.2 #burn rate coefficient (mm/s) (this can actually vary with pressure so if there is data for that, this can be included)
    n = 0.688 #burn rate exponent (unitless) (this can actually vary with pressure so if there is data for that, this can be included)
    kv = 0 # Propellant erosive burning velocity coefficient (unitless)
    ci = 1 #not sure
    osi = 0 #not sure
    ei = 1 #not sure
    rhopgrain = 1.78505 #density of propellant (kg/m^3)
    dt0 = 16.79248 #inital nozzle throat diameter (mm)
    nc = 0.95 #combustion effiecency. Use 0.95 for fine grain and 0.93 for coarse grain. (unitless)
    temp0 = 1710 #ideal combustion temperature (kelvin)
    k = 1.1308 #ratio of specific heats, mixture (unitless)
    m = 42.42 #effective moleculer weight of products (kg/kmol)
    d1 = 20 #grain core diamter (mm)
    d2 = 69 #grain outside diamter (mm)
    l = 460 #grain total length (mm)

    #changing variables
    xi = 0 #Web regression (mm)
    xinc = 0.0293764988009578 #x increment (mm)

    
    #derivable variables
    tweb0 = (d2 - d1)/2 #initial web thickness (mm)
    vc = (math.pi / 4) * (dc**2) * lc #volume of chamber (mm^3)
    mgrain = [rhopgrain*(((math.pi/4) * (d2**2 -d1**2)*l)/(1000**3))/(1000**2)]
    rat = 8314 / m #specific gas constant (J/kg-K)
    p0 = [patm] #absoulte chamber pressure (this starts as local atmospheric pressure) (MPa) 
    t = [0, xinc] #time since start of burn (s)
    masssto = [0]


    first = True
    #pressure calc
    while p0[-1]-patm > 0 or first: #run through burn until gauge pressure is 0 (burn has ended)
        first = False
        #iteration updates
        xi = xi + xinc
        d1 = d1 +(ci * 2 *xinc)
        d2 = d2 -(osi *2 * xinc)
        l = l-(ei*n*2*xinc)
        tweb = (d2-d1)/2 #web thickness (mm)
        at = (math.pi/4) * (dt0+(((erate)*(tweb0-tweb))/tweb0))**2 #Throat area (mm^2) stays constant without an erosion rate
        astar = at/(1000**2) #Nozzle critical passage area (m^2)
        aduct = ((math.pi/4) * (dc**2)) - (math.pi/4)*(d2**2 - d1**2) #difference in chamber and grain cross-sectional area (mm^2)
       
        g = 0 #erosive burning factor (unitless)
        if gstar - (aduct/at) > 0:
            g = gstar - (aduct / at)
       
        r = (1 + (kv*g)) * a * (p0[-1]**n) #propellant burn rate (mm/s)
        vgrain = ((math.pi/4) * (d2**2 -d1**2)*l)/(1000**3) #grain volume (m^3)
        vfree = (vc / 1000**3) - vgrain #free volume in chamber (m^3)
        mgrain.append((rhopgrain*vgrain)/(1000**2)) #mass of propellent remaining (kg)
        mgen = (mgrain[-2] - mgrain[-1] ) / (t[-1] - t[-2] ) #mass generation rate of combustion products (kg/s)
        mnoz = (p0[-1] - patm)*1000000*astar/math.sqrt(temp0*rat*nc)*math.sqrt(k)*((2/(k+1))**((k+1)/2/(k-1))) #mass flow through nozzle (kg/s)

        if mgen < mnoz:
            if p0[-1] < 0:
             mnoz = 0

        msto = mgen - mnoz #Mass storage rate of combustion products (in chamber) (kg/s)

        if msto*(t[-1] - t[-2]) + masssto[-1] < 0:
            masssto.append(0)
        else:
            masssto.append(msto*(t[-1] - t[-2]) + masssto[-1])

        rhoprod = masssto[-1] / vfree #density of cobustion products in chamber (kg/m^3)

        p0.append(((rhoprod * rat * temp0*nc) + patm*1000000)/1000000) #absolute chamber pressure (MPa)

        
        t.append(xinc/(r)+t[-1])

    return [t,p0]

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

