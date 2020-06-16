# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:01:15 2020

This program convets a framework cif file into a vasp readible POSCAR and
generates "N" random orientations of a molecule (within POSCAR format) in the 
vacuum gap between framework layers.

This is the first program for the 'materials_bandgap_workflow' project.

@author: chr218
"""

import sys
import numpy as np
import openbabel
from ase.io  import vasp, cif
from ase import Atoms
from scipy import constants,optimize
from ase.visualize import view

'''
constants
'''



'''
subroutines & functions
'''

    

'''
'get_framework' converts the framework cif file into vasp readible CONTCAR and
returns an ASE type atoms object.

ASE sucks at converting '.cif' to 'CONTCAR'. Openbabel is much better.

This function first converts the '.cif' into an openbabel type object, 
saves as 'CONTCAR' type, then reads this 'CONTCAR' using ASE. 

Highly inefficient! But fix later.

Finally, it returns a repeated atoms object defining the framework.
'''
def get_framework(fname,rep):
    #Define openbabel file conversion type. "from" and "to"
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("cif", "CONTCAR")
    
    #Define openbabel molecule type object.
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, fname)
    
    #Convert and write
    obConversion.WriteFile(mol, 'CONTCAR_frame')

    #Read using ASE, generate atoms type object.
    frmwrk_atoms = vasp.read_vasp('CONTCAR_frame')
    
    #Repeat atoms.
    if len(rep) > 0:
        frmwrk_rep = frmwrk_atoms.repeat(rep)
    
    
    return frmwrk_rep


'''
'gen_orient' generates randomly oriented molecules within the framework of
interest.
'''    
def gen_orient(mol,frm,N_images):
    i = 0
    while i < N_images:
        #np.random.seed(i)
        
        #Create copy of atoms object 'mol'.
        mol_cp = mol.copy()

        #Randomly generate rotation about x,y,z axis for gas mol placement.
        rot_gen = np.multiply(np.random.random_sample([3]),np.array([180,180,180]))
        mol_cp.euler_rotate(phi=rot_gen[0],theta=rot_gen[1],psi=rot_gen[2],center=(0,0,0))
        
        #Randomly generate center of mass position for gas mol to be placed.

        com_gen = np.random.random_sample([3])

        #Replicate com_gen vector along axis= 0 to form array.
        com_gen = np.tile(com_gen,(len(mol.get_chemical_symbols()),1))
        
        #Translate gas mol to new com position and wrap atoms outstanding from unit cell.
        mol_cp.set_scaled_positions(np.add(mol.get_scaled_positions(),com_gen))
        mol_cp.wrap()
    
        #Create system 'frm_mol' by appending mol atoms object to framework atoms object.
        mol_frm = mol_cp+frm
        
        #Ensure new com position does not overlap with framework atoms.
        atoms_far_enough_criteria = 0
        for j,atom in enumerate(mol_frm[0:mol.get_global_number_of_atoms()]):
            dist = mol_frm.get_distances(j,range(mol.get_global_number_of_atoms(),mol_frm.get_global_number_of_atoms()))
            if np.amin(dist) > 5.0:#Adjust this for molecule placement. Helps avoid overlap!
                atoms_far_enough_criteria += 1
                
        if atoms_far_enough_criteria == mol.get_global_number_of_atoms():
            f_out = ('CONTCAR_%s' % str(i))
            vasp.write_vasp(f_out,mol_cp+frm)
            view(mol_frm)
            i += 1
            continue

'''
main program begins here
'''
#Generate atoms objects for molecule and framework.
gas_mol = vasp.read_vasp('methane\CONTCAR')
frmwrk = get_framework('frameworks\MoS2_mp-1027525_computed.cif',(3,4,1))

#Set molecule's unit cell to be that of framework's unit cell.
gas_mol.set_cell(frmwrk.get_cell_lengths_and_angles())

#Center the molecule's center of mass at the origin of the framework's unit cell.
gas_mol.center(about=0)

#Generate random orientations wihtin framework.
gen_orient(gas_mol,frmwrk,5)



