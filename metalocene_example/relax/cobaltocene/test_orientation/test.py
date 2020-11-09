import sys
import numpy as np
import openbabel
from ase.io  import vasp, cif
from ase import Atoms
from scipy import constants,optimize
from ase.visualize import view

gas_mol = vasp.read_vasp('POSCAR_unwrapped')
gas_mol.center(about=0)

view(gas_mol)

mol_cp = gas_mol.copy()

rot_gen = np.multiply(np.random.random_sample([3]),np.array([180,180,180]))



mol_cp.euler_rotate(phi=rot_gen[0],theta=rot_gen[1],psi=rot_gen[2],center=('COM'))



view(mol_cp)
