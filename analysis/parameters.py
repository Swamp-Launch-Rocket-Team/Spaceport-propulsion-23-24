import math

R_univ = 8.314

class Params:
    class nozzle:
        def __init__(self):
            self.init = False
            self.throat_diam = 0
            self.throat_area = 0
        
    class prop:
        def __init__(self):
            self.init = False
            self.burn_rate = 0
            self.burn_temp = 0
            self.dens = 0
            self.char_velo = 0
            self.heat_ratio = 0
            self.molar_mass = 0
            self.burn_coef = 0
            self.burn_exp = 0
            self.KN = 0
            self.Di = 0
            self.Do = 0
            self.L = 0
            self.R = 0

            

    def __init__(self):
        self.nozzle = Params.nozzle()
        self.prop = Params.prop()
        self.dt = 0
 
        