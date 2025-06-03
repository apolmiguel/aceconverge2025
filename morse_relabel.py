import os 
import sys
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from y2ace_funcs import *
import itertools

## FUNCTIONS ## 
def morse_potential(r, a=8.57696795, b=2.00550983, c=1.28317311, d=-8.57696795):
    return a * (1 - np.exp(-b * (r - c)))**2 + d

def force_morse(r, a=8.57696795, b=2.00550983, c=1.28317311):
    """ Force kernel from the morse potential """
    return 2 * a * b * (1 - np.exp(-b * (r - c))) * np.exp(-b * (r - c))

# From PANNA # 
def make_replicas(positions, lattice_vectors, Rc, pbc=[True, True, True]):


    """ Function computing directly all the replicas needed to contain
        all atom images within the cutoff.
        Also wraps original positions for convenience

    Parameters
    ----------
    lattice_vectors: lattice vectors as matrix (a1, a2, a3)
    Rc: maximum radial cutoff that you want to take in to account
    pbc: Normally pbc are recovered from lattice vector,
         if in the lattice_vectors a direction is set to zero
         then no pbc is applied in that direction.
         This argument allow you to turn off specific directions
         in the case where that specific direction has a lattice vector
         greater then zero.
         To achieve this pass an array of 3 logical value (one for each
         direction). False value turn off that specific direction.
         Default is true in every direction.
         eg. pbc = [True, False, True] => pbc along a1 and a3

    Returns
    -------
    positions: the wrapped atomic positions
    replicas: numpy array of displacement vectors for each replica.
    """
    lattice_vectors = np.asarray(lattice_vectors)
    # wrap the atoms in the unit cell
    if (lattice_vectors == 0).all():
        pass
    else:
        crystal_coord = positions @ np.linalg.inv(lattice_vectors)
        crystal_coord = crystal_coord % 1.0
        positions = crystal_coord @ lattice_vectors

    # compute how many cell replica we need in each direction
    max_indices = replicas_max_idx(lattice_vectors, Rc, pbc)
    l_max, m_max, n_max = max_indices
    l_list = range(-l_max, l_max + 1)
    m_list = range(-m_max, m_max + 1)
    n_list = range(-n_max, n_max + 1)
    # create a matrix with all the idx of the extra cell we need
    replicas = np.asarray(list(itertools.product(l_list, m_list, n_list)))

    # translation vectors needed to map each atom in the unit cell to its replica
    # @ means matrix multiplication
    replicas = replicas @ lattice_vectors

    # Keep only the replicas that are actually inside cutoff
    # Since we use all atoms in the tf code, this can be a good speedup
    corners = np.asarray(list(itertools.product([0,1],[0,1],[0,1]))) @ lattice_vectors
    keep_replicas = []
    for r in replicas:
        distances = np.sum((corners[:,np.newaxis,:] - \
                            corners[np.newaxis,:,:] - \
                            r[np.newaxis,np.newaxis,:])**2,axis=2)
        if np.min(distances)<Rc**2:
            keep_replicas.append(r)

    return positions, np.asarray(keep_replicas)

## DATA LOAD ## 
filedir = sys.argv[1]
dset = load_gzip(filedir)
rcut = float(sys.argv[2])

for ind,struc in enumerate(dset):

    lvecs = dset['ase_atoms'][ind].get_cell()
    pos = dset['ase_atoms'][ind].get_positions()
    nats = len(pos)
    print('Number of atoms:',nats)
    print('Corrected energy:',dset['energy_corrected'][ind])
    
    _, replicas = make_replicas(pos, lvecs, rcut)
    posi = np.reshape(pos,(nats,1,1,3))
    posj = np.reshape(pos,(1,nats,1,3))
    reps = np.reshape(replicas,(1,1,-1,3))
    rij_vec = posj + reps - posi
    rij = np.linalg.norm(rij_vec, axis=3)
    rij = np.where((rij > rcut) | (rij < 1e-8), np.inf, rij)
    print('Pre-corrected energy:', dset['energy'][ind])
    print('Morse potential energy:',0.5*np.sum(morse_potential(rij)))
    energy = 0.5 * np.sum(morse_potential(rij))
    forces_struc = (1/rij * force_morse(rij))[:,:,:,np.newaxis] * rij_vec
    forces = np.sum(forces_struc, axis=(1,2))

    # Replacing energies and forces in the dataset
    dset['energy'].iloc[ind] = energy
    dset['energy_corrected'].iloc[ind] = energy
    dset['forces'].iloc[ind] = forces.tolist()

## REBUILD DATASET ## 
prefix = filedir[:-10]
dset.to_pickle(prefix + '_morse.pckl.gzip', compression='gzip', protocol=4)
print('Dataset saved as:', prefix + '_morse.pckl.gzip')
    


