import numpy as np
import pandas as p
import sys 

from y2ace_funcs import *

# Variables
property_string = 'Properties=species:S:1:pos:R:3:force:R:3:tags:I:1' 
cell_string = 'unit_cell=conventional'
config_string = 'config_type=shaiducarbon'
gzip_file = sys.argv[1]
outname = sys.argv[2]

def gzip2extxyz(gzipfile, outname, verbose=False):
    """
    Convert a gzip file containing atomic positions and forces to an extxyz file.
    
    Parameters:
    gzipfile (str): Path to the input gzip file.
    outdir (str): Directory where the output extxyz file will be saved.
    outname (str): Name of the output extxyz file (without extension).
    
    Returns:
    None
    """

    dataset = load_gzip(gzipfile)

    with open(outname, 'w') as file:
        for i, entry in enumerate(dataset['ase_atoms']):
            natom = str(len(entry))
            lattice_string = 'Lattice=' + '\"' + " ".join(f"{j}" for j in entry.get_cell()[:].flatten().tolist()) + '\"'
            energy_string = 'energy=' + str(dataset['energy'].iloc[i])
            pbc_string = 'pbc=' + '\"' + " ".join(f"{j}" for j in np.where(entry.get_pbc(),'T','F')) + '\"'
            file.write(f"{natom}\n{lattice_string} {property_string} {config_string} {energy_string} {cell_string} {pbc_string}\n")
            if verbose:
                print(natom,'\n' + lattice_string + ' ' + property_string + ' ' + config_string + ' ' + energy_string + ' ' + cell_string + ' ' + pbc_string)

            species = [[k] for k in entry.get_chemical_symbols()]
            pos = entry.get_positions().tolist()
            forces = dataset['forces'].iloc[i]
            for m,coords in enumerate(pos):
                for n,alpha in enumerate(coords):
                    pos[m][n] = f"{float(pos[m][n]):.8f}"
                    forces[m][n] = f"{float(forces[m][n]):.8f}"
            per_atom = [sp + p + f + ['0'] for sp, p, f in zip(species, pos, forces)]
            for entry in per_atom:
                file.write('\t\t'.join(f"{i}" for i in entry) + '\n')
                if verbose:
                    print('\t'.join(f"{k}" for k in entry))

if __name__ == "__main__":
    gzip2extxyz(gzip_file, outname, verbose=False)
    print('Conversion complete. Output saved to', outname)