import os
import numpy as np 
from ase.io import read, write
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# from ase import Atoms

# from y2ace_funcs import *

## Data load ## 
dataset = read('datasets/Si2.xyz',index=':') # ':' gets all frames
print(len(dataset))
print(dataset[0])
e1b = -7881.32677981122 # isolated energy of Si atom (chong2025)

## Si2 dimer reconstruction ## 
elist = []
distlist = []
for struc in dataset[:]:
    positions = struc.get_positions()
    energy = struc.get_potential_energy()
    elist.append(energy)
    distlist.append(np.linalg.norm(positions[0] - positions[1]))
elist = np.array(elist) - e1b * 2
distlist = np.array(distlist)
data = np.column_stack((distlist, elist))
print(data)
np.savetxt('datasets/Si2_energy.dat', data, fmt='%.8f')


## Extract positions, generic ##
# distlist = []
# for struc in dataset[:]:
#     positions = struc.get_positions() 
#     posi = positions[:, np.newaxis, :]
#     posj = positions[np.newaxis, :, :]
#     rij_vec = posi - posj
#     rij = np.linalg.norm(rij_vec, axis=-1)
#     distlist.append(rij[np.where(rij > 1e-8)])
# distlist = np.array(distlist) # actual count is halved 
# print(np.histogram(distlist.flatten()))
    