import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
from ase import Atoms
from pyace import PyACECalculator
import matplotlib.cm as cm

import os
import json
import numpy as np
import pandas as pd
from ase import Atoms
from pyace import PyACECalculator
from scipy.optimize import curve_fit
from pbc import orthogonal_vector,replicas_max_idx,make_replicas

### FUNCTIONS ### 

def en_ase(r,potential_dir,potential_name):
    '''
    Returns:
    -------
    Total dimer/two-body energy array as a function of distance r.
    Potential is called by PyACECalculator.

    Inputs:
	-------
    1. r - (1,ndata) numpy array of positions. ndata is the length of the array input 
    2. potential_dir - Directory where the .yaml files of the potential (output or interim) can be found. 
    3. potential_name - Filename of chosen potential.
    '''
    ndata = len(r)
    en_array = np.empty([ndata,1])
    calc = PyACECalculator(potential_dir + potential_name)
    for i,rval in enumerate(r):
        atoms = Atoms('C2', positions=[[0,0,0],[rval,0,0]])
        atoms.set_calculator(calc)
        en_array[i] = atoms.get_potential_energy()
    return en_array

def collect_from_json(filedict):
    '''
    Returns:
    --------    
    spec_label: array of element labels per atom; new in y2ace  
    config_energy: array of total energy of configuration and units
    lattvecs: lattice vectors
    positions: position vectors of each atom
    forces: force vectors of each atom
    
    Parameters:
    -----------
    filedict: dictionary loaded from .example file

    '''
    length_units = filedict['unit_of_length']
    config_energy = filedict['energy']
    lattvecs = np.array(filedict['lattice_vectors'])
    species_data = filedict['atoms']
    spec_label = [sval[1] for sval in species_data] # species label list
    pos_forces = [] # parsing position, and forces  
    for i,sval in enumerate(species_data):
        pos_forces.append(sval[2:])
    pos_forces = np.array(pos_forces) # pos_forces[atom][pos,force][x,y,z]
	
    positions,forces = pos_forces[:,0],pos_forces[:,1]

    return spec_label,config_energy,length_units,lattvecs,positions,forces

def load_gzip(gzip_file):
	'''
	Input/s:
	--------
	gzip_file: path to .gzip file, including the filename.

	Returns:
	--------
	Pandas dataframe of gzip_file

	'''
	return pd.read_pickle(gzip_file,compression='gzip')

def energy_per_atom(data_dir):
	'''
	Returns: 
	--------
	Numpy array of energy per atom for a given dataframe.

	Input:
	------
	data_dir: String directory of the pacemaker.pckl.gzip file.
	'''
	data = pd.read_pickle(data_dir,compression='gzip')
	headers = list(data)
	
	# print error if both energy and ase_atoms not in headers
	if 'energy' not in headers:
	    print('Error: energy not in headers')
	    return
	if 'ase_atoms' not in headers:
	    print('Error: ase_atoms not in headers')
	    return

	energy_array = data['energy']
	natom_array = np.empty(len(data['ase_atoms']))
	for ind,conf in enumerate(data['ase_atoms']):
	 natom_array[ind] = len(conf.get_positions())

	energy_per_atom = energy_array/natom_array

	return energy_per_atom

def distmat_builder(array1,array2):
    ''' 
    Returns:
    --------
    Generalized distance matrix with dimensions (len(array1),len(array2),3)

    Input/s:
    -----------
    array1,array2: numpy arrays of generic dimensions. Must be broadcastable. 

    '''
    return np.sqrt(np.sum((array1-array2)**2,axis=-1))


def strucs2dists(ase_atoms,rcut,warn=True):
    '''
    Returns:
    --------
    1D numpy array of distances between atoms within a cutoff in a given ase.Atoms structure.
    ! Caution !: Distances are counted twice. Keep in mind when using.
    
    Input/s:
    --------
    ase_atoms - ase.Atoms structure where interatomic distances are measured. 
    rcut - Cutoff distance. Distances r not counted beyond r > rcut. 

    Prerequisites:
    --------------
    1. ase package (https://wiki.fysik.dtu.dk/ase/) imported to be able to read the ase.Atoms object. 
    2. Environment where this is run should have access to PANNA's (https://pannadevs.gitlab.io/pannadoc/) gvector/pbc.py to make the replicas.

    NB: Function loosely based off distmat_builder() and remove_e2b() from twobody_remover.py of https://github.com/apolmiguel/sissa_y1proj.
    '''

    ### Get position and lattice vector arrays from ase.Atoms structure ### 
    positions = ase_atoms.get_positions()    
    lattvecs = ase_atoms.get_cell()[:3] 

    ### Generate array of replicas within rcut using PANNA's pbc.py ### 
    atoms,reps = make_replicas(positions,lattvecs,Rc=rcut)
    images = atoms[:,np.newaxis,:] + reps[np.newaxis,:,:]

    ### Build distance matrix between all atoms and images. Filter self-distances and those outside of the cutoff. ###
    distmat = distmat_builder(atoms[:,np.newaxis,np.newaxis,:],images[np.newaxis,:,:,:])
    dist_filt = distmat[(distmat!=0) & (distmat<=rcut)]


    ### Output vector and include warning about the distances being counted twice. ### 
    if warn:
      print('Reminder: strucs2dists() output will have the distances counted twice.')
    else:
      pass
    return dist_filt

def get_mindist(gzip_file,rcut=6):
    '''
    Returns:
    --------
    Float value of minimum distance between atoms in a given ase.Atoms structure.
    
    Input/s:
    --------
    gzip_file: path to .gzip file, including the filename.
    rcut: Cutoff distance. Distances r not counted beyond r > rcut.
    '''
    gzip = load_gzip(gzip_file)
    mindists = []
    for structure in gzip['ase_atoms']:
	    dists = strucs2dists(structure,rcut,warn=False)
	    mindists.append(np.min(dists))
    return np.min(mindists)


### From monkeybars.ipynb (tester notebook) ###

def load_rmse(filedir,trainflag=True,mevflag=True):
    '''
    Input/s:
    --------
    filedir - Directory where the "ladder_rmse.dat" file is.
    trainflag
        True: obtains 
    mevflag - Whether the energy and forces get converted from eV to meV.

    Returns:
    --------
    Numpy array of loaded rmse file.
    '''
    if trainflag==True: 
        rmse_list = np.loadtxt(filedir+'/train_rmse.dat').T
    else:
        rmse_list = np.loadtxt(filedir+'/test_rmse.dat').T
    if mevflag==True:
        rmse_list[1] *= 1000 # eV to meV
        rmse_list[2] *= 1000 # eV/A to meV/A
    return rmse_list

## Updated version to include the potential and remove flag variables. ## 
def load_error_dimerpot(filedir,dist_array,potname='output_potential.yaml'):
    '''
    Returns:
    --------
    (1) Numpy arrays of training (train_rmse.dat) and validation (test_rmse.dat) files.
    (2) Numpy array of potentials evaluated at dist_array distances. 
    
    Input/s:
    --------
    filedir - Directory where the error and potential files are.
    dist_array - Numpy array of dimer distances used as arguments for the potname.yaml ASE potential.
    potname - String prefix of .yaml file to extract. 'output_potential' by default, but can get interim potentials.
    '''
    train_rmse = np.loadtxt(filedir+'train_rmse.dat').T
    test_rmse = np.loadtxt(filedir+'test_rmse.dat').T
    ase_pot = en_ase(dist_array,filedir,potname)
    return train_rmse,test_rmse,ase_pot

## Function that loads load_error_dimerpot output to lists iteratively ## 
def multiload_error_dimerpot(filedir_main,dist_array,potname='output_potential.yaml',body_range=np.arange(2,5)):
    '''
    Returns:
    --------
    Lists of load_error_dimerpot output according to each value of body_range.

    Input/s:
    --------
    Almost identical to load_error_dimerpot except for
    filedir_main - Directory where the other body order (ex. b_order3) folders are. b_order+str(i) should contain the files for the errors and potentials.
    body_range - integer list of which body order (default - 2, 3, 4) to include.
    Note: file directory of the files should be in the format of filedir+'border'+str(body_range[ind]).
    '''
    trainerrorlist = []
    testerrorlist = []
    asepotlist = []
    for bd_ind in body_range:
        trainerror,testerror,asepot = load_error_dimerpot(filedir_main+'b_order'+str(bd_ind)+'/',dist_array,potname)
        trainerrorlist.append(trainerror)
        testerrorlist.append(testerror)
        asepotlist.append(asepot)
    return trainerrorlist,testerrorlist,asepotlist

## Function that extracts elements in the training and validation error lists generated by multiload_error_dimerpot per body order ## 
def extract_besterror(trainerrorlist,testerrorlist,bodyorders,whicherror=3):
    '''
    Returns:
    -------- 
    train/testerrors -  Numpy list of best training/validation errors per body order.
    
    Input/s:
    --------
    train/testerrorlist - Numpy arrays of test/validation errors made by multiload_error_dimerpot.
    bodyorders - Integer value of how many body orders were considered.
    whicherror - Integer indicator of which metric to get: 0 - nfuncs (not important) , 1 - RMSE energy, 2 - RMSE force comp, 3 (or -1) - loss function.
    '''
    trainerrors = [] # empty list containers 
    testerrors = []
    for i in np.arange(bodyorders):
        trainerrors.append(trainerrorlist[i][whicherror][-1]) 
        testerrors.append(testerrorlist[i][whicherror][-1])
    return np.array(trainerrors),np.array(testerrors)

## Making colormap values match with alpha ## 
def cm_w_alpha(hexcode,Npts,bias=0.1,slope=0.9):
    '''
    Returns:
    --------
    cm_w_a - A [(hexcode,alpha)] list of tuples of len(Npts) for scatter plotting a single color generated from colormaps.cm with variable alphas.
             Can serve as a substitute for a cm.var specific on Npts.
    
    Inputs:
    -------
    # cmapvar - a matplotlib.cm.colormap variable (array of len 4 [r,g,b,alpha] representing a single color.
    hexcode - a hexcode string of the single color I wish to enter.
    Npts - number of points for the desired plot.
    bias - Float between 0 and 1 indicating the starting alpha/transparency value.
    slope - Float between 0 and 1 - bias indicating how sharp the alpha/transparency changes.
    '''

    facecolors = [hexcode]*Npts
    facealphas = [bias+slope*i/Npts for i in range(Npts)]
    cm_w_a = list(zip(facecolors,facealphas))
    return cm_w_a


### MAIN ###
if __name__ == '__main__':
	print(get_mindist('../datasets/Carbon_full/Tr50k_n.pckl.gzip',rcut=6))



### Plotting functions needed for building the files here.### 

## Learning rate ##
def plot_learnrates(datasets,trainrates,valrates,plot_title,filename,val_flag=True,cmap_flag=False):
	'''
	Returns:
	--------
	Training and validation learn ratess of energies, force components, and loss of the datasets.

	Parameters:
	-----------
	datasets: String list of dataset labels. Len should be consistent with trainrates and valrates.
	trainrates: List of numpy arrays of training learning rates (following format of lrnrate.dat).
	valrates: List of numpy arrays of validation learning rates (following format of vallrn.dat).
	plot_title: String title for the plot.
	filename: String filename for the plot.
	val_flag: Flag whether validation rates are plotted (True). Default true.
    cmap_flag: Flag whether perceptually-uniform color map is used (True). Default false.
	'''

	fig1 = plt.figure(figsize=(14,7),dpi=200)
	fig1.subplots_adjust(top=0.95,wspace=0.2,hspace=0.2)
	fig1.suptitle(plot_title)
	axle = fig1.add_subplot(221); axlf = fig1.add_subplot(223); axl = fig1.add_subplot(122)
	axset = [axl,axle,axlf]
	axylabels = ['Optimized loss','RMSE_e (meV/at)', 'RMSE_fcomp (meV/A)']


	if cmap_flag:
		clist = [cm.plasma(val) for val in np.linspace(0,0.9,len(datasets))[::-1]]
	else:
		clist = ['C'+str(i) for i in range(len(datasets))]

	for ind,ax in enumerate(axset):
		if val_flag:
			for jnd,dsl in enumerate(datasets):
				ax.plot(trainrates[jnd][0],trainrates[jnd][ind+1],lw=2,linestyle='--',color=clist[jnd])
				ax.plot(valrates[jnd][0],valrates[jnd][ind+1],lw=2,label=dsl,color=clist[jnd])
		else:
			for jnd,dsl in enumerate(datasets):
				ax.plot(trainrates[jnd][0],trainrates[jnd][ind+1],lw=2,label=dsl,linestyle='--',color=clist[jnd])
		# ax.grid(axis='y',which='both')
		ax.set_xscale('log'); ax.set_yscale('log')
		ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
		#ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
		ax.set_ylabel(axylabels[ind])
		ax.set_xlabel('Iterations')
		ax.legend()
	fig1.savefig(filename,bbox_inches='tight')
	print('Plot saved:',filename)

def plot_dimercurves(potdirs,potname,datasets,plot_title,filename,mindist,rcut,ymin=-10,ymax=5,cmap_flag=False):
	'''
	Returns:
	--------
	Plot of dimer curve energies of the potentials (potname) in each of the directories in potdirs as a function of input distance vector, r.

	Parameters:
	-----------
	r	   : Numpy array of input distances.
	potdirs: String list of directories containing the runs.
	potname: String name of the potential; default is interim_potential_0.yaml (for prematurely-stopped runs).
	mindist: Minimum distance found in the datasets to be plotted (set for xlim).
	rcut   : Cutoff distance for the potential.
    ylim   : Tuple that dictates the y-limits of the plot.
	datasets: String list of dataset labels.
	plot_title: String title for the plot.
	filename: String filename for the saved figure.
    cmap_flag: Flag whether perceptually-uniform color map is used (True). Default false.

	Prerequisites:
	--------------
	- Dimer curve data from QE calculations in datafiles/es01_dimer.dat.
	- PyACECalculator from pyace.
	- Atoms from ASE.
	- en_ase function from y2ace_funcs.
	'''
	if cmap_flag:
		clist = [cm.viridis(val) for val in np.linspace(0,0.7,len(datasets))]
	else:
		clist = ['C'+str(i) for i in range(len(datasets))]
	# Data: Distances and dimer energies from QE #
	r = np.linspace(1,rcut+1,200)
	en_qe = np.loadtxt('datafiles/es01_dimer250.dat').T
	en_qe[0] *= 0.529177 # bohr to Angstrom
	en_qe[1] -= -18.03977639 * 2 # Subtracting the energy of two isolated atoms.
	en_qe[1] *= 13.605703976 # Ry to eV

	# Plotting #
	g,gax = plt.subplots(1,1,figsize=(7,3.5),dpi=200)
	# g.suptitle(plot_title)
	g.subplots_adjust(top=0.95)
	# gax.scatter(*en_qe,s=1,color='k',label='DFT')
	gax.scatter(*en_qe,s=1,color='k')
	if cmap_flag:
		for i,dir in enumerate(potdirs):
			gax.plot(r,en_ase(r,dir,potname),label=datasets[i],color=clist[i])
	else:
		for i,dir in enumerate(potdirs):
			gax.plot(r,en_ase(r,dir,potname),label=datasets[i])        
	gax.set_xlabel('Distance (A)')
	gax.set_ylabel('Energy (eV)')
	gax.set_xlim(mindist,rcut+1)
	gax.set_ylim(ymin,ymax)
	gax.legend(fontsize=8)
	g.savefig(filename,bbox_inches='tight')
	print('Plot saved:',filename)
	

