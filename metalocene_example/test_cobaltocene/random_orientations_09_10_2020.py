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
import subprocess
#from openbabel import openbabel
from ase.io  import vasp, cif
from ase import Atoms
from scipy import constants,optimize
#from ase.visualize import view

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
    #obConversion = openbabel.OBConversion()
    #obConversion.SetInAndOutFormats("cif", "CONTCAR")
    
    #Define openbabel molecule type object.
    #mol = openbabel.OBMol()
    #obConversion.ReadFile(mol, fname)
    
    #Convert and write
    #obConversion.WriteFile(mol, 'CONTCAR_frame')

    #Read using ASE, generate atoms type object.
    frmwrk_atoms = vasp.read_vasp('/home/chr218/materials_bandgap_workflow/materials_bandgap_workflow-master/example/frameworks/adjusted_vaccum/CONTCAR_frame')
    
    #Repeat atoms.
    if len(rep) > 0:
        frmwrk_rep = frmwrk_atoms.repeat(rep)
    
    
    return frmwrk_rep


#Generate POTCAR for CONTCAR
def make_POTCAR(atoms):
        symbol_list = []
        for atom in atoms:
            symbol_list.append(atom.symbol)

        symbol_arr = np.asarray(symbol_list)
        POTCAR_symbols = np.unique(symbol_arr)
        #POTCAR_symbols = (symbol_arr)

        print('Generating POTCAR for atoms: \n')
        for sym in POTCAR_symbols:
                print(sym)
                subprocess.call(['cat /home/chr218/bin/POTCAR_files/PBE_52/POTCAR_'+sym+' >> POTCAR'],shell=True)

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

        #Randomly generate rotation about COM axis for gas mol placement.
        rot_gen = np.multiply(np.random.random_sample([3]),np.array([180,180,180]))

        mol_cp.euler_rotate(phi=rot_gen[0],theta=rot_gen[1],psi=rot_gen[2],center=('COM'))
        
        #Randomly generate center of mass position for gas mol to be placed.
        com_gen = np.random.random_sample([3])

        #Replicate com_gen vector along axis= 0 to form array.
        com_gen = np.tile(com_gen,(len(mol_cp.get_chemical_symbols()),1))
        
        #Translate gas mol to new com position and wrap atoms outstanding from unit cell.
        mol_cp.set_scaled_positions(np.add(mol_cp.get_scaled_positions(),com_gen))
        mol_cp.wrap()
    
        #Create system 'frm_mol' by appending mol atoms object to framework atoms object.
        mol_frm = mol_cp+frm
        
        
        #Ensure new com position does not overlap with framework atoms.
        atoms_far_enough_criteria = 0
        for j,atom in enumerate(mol_frm[0:mol.get_global_number_of_atoms()]):
            dist = mol_frm.get_distances(j,range(mol.get_global_number_of_atoms(),mol_frm.get_global_number_of_atoms()))
            if np.amin(dist) > 1.0:#Adjust this for molecule placement. Helps avoid overlap!
                atoms_far_enough_criteria += 1
                
        if atoms_far_enough_criteria == mol.get_global_number_of_atoms():
            f_out_dir = ('CONTCAR_%s' % str(i))
            f_out = ('POSCAR')
            vasp.write_vasp(f_out,(mol_cp+frm),sort=True)
            #view(mol_frm) #uncomment to view formed CONTCARs within ASE gui.i

            #Generate POTCAR file.
            make_POTCAR(mol_cp+frm)

            #Generate INCAR file. (Just pull it from bin for now- FIX LATER)
            subprocess.call("cp ~/bin/INCARS_other_functionals/DFT-D2/INCAR .",shell=True)
            #Generate submission script, again just pull from bin for now.
            subprocess.call("cp ~/bin/submit_vasp_gam.sh .",shell=True)
            #Generate KPOINTS file, again just pull from bin for now.
            subprocess.call("cp ~/bin/KPOINTS .",shell=True)

            subprocess.call("mkdir "+f_out_dir,shell=True)
            subprocess.call("mv KPOINTS POTCAR submit_vasp_gam.sh INCAR POSCAR "+f_out_dir,shell=True)
            i += 1
            continue

'''
main program begins here
'''
#Generate atoms objects for molecule and framework.
gas_mol = vasp.read_vasp("./POSCAR_unwrapped")
frmwrk = get_framework("/home/chr218/materials_bandgap_workflow/materials_bandgap_workflow-master/example/frameworks/MoS2_mp-1023939_computed.cif",(3,4,1))

#Set molecule's unit cell to be that of framework's unit cell.
gas_mol.set_cell(frmwrk.get_cell_lengths_and_angles())

#Center the molecule's center of mass at the origin of the framework's unit cell.
gas_mol.center(about=0)

#Generate random orientations wihtin framework.
gen_orient(gas_mol,frmwrk,10)



