# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pyswmm import Simulation, Nodes, Links

fname = 'C:/Users/ga.riano949/Documents/GitHub/mpc/Paper ACC 2018/models/parameter_estimation/test pyswmm/4states_simple.inp'

f = []
with Simulation(fname) as sim:
    links = Links(sim)
    nodes = Nodes(sim)
    for step in sim:
        ff = links['T1'].flow
        print(ff)
        f.append(ff)
        
        
    
plt.plot(f)
plt.show()