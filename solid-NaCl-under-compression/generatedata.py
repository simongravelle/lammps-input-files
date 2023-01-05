#!/usr/bin/env python
# coding: utf-8

# # Create a data file for LAMMPS
# ## The system consist of a NaCl crystal wall in contact with a gas of Ethanol molecules
# ## Import libraries
import numpy as np
import copy

# ## parameter choice
dnacl = 2.84 # Na-Cl typical distance
dx, dy, dz = dnacl, dnacl, dnacl
nx, ny, nz = 10, 10, 10
Lx, Ly, Lz = nx*dx, ny*dy, nz*dz
print('The desired Na-Cl layer dimensions are '+str(Lx)+' A x ' +str(Ly)+' A x '+str(Lz)+' A')


# ## Lammps box size
txlo, txhi = -Lx/2, Lx/2
tylo, tyhi = -Ly/2, Ly/2
tzlo, tzhi = -Lz/2, Lz/2

# ## create the NaCl wall
basestructure = np.loadtxt('../../shared/NaCl/Position.dat') # import the 8 basic atoms to replicates
# replicate the initial structure
naclwall = copy.deepcopy(basestructure)
for xx in np.arange(txlo+dx/2,txhi,2*dx):
    for yy in np.arange(tylo+dy/2,tyhi,2*dy):
        for zz in np.arange(tzlo+dz/2,tzhi,2*dz):
            naclwall = np.append(naclwall,basestructure+[0,0,0,0,xx,yy,zz], axis=0)
naclwall = naclwall[8:]
# renumber atoms ids
for n in range(len(naclwall)):
    naclwall[n,0] = np.int64(n+1)

# ## create data lammps file
cptatm = 0
cptmol = 1
atoms = np.zeros((1000000,7))
# ## place the NaCl
for m in naclwall:
    atoms[cptatm] = m[0], cptmol, m[2], m[3], m[4], m[5], m[6]
    cptatm += 1

# ## remove excess lines
atoms = atoms[0:cptatm]       

# ## write LAMMPS data file
f = open("data.lammps", "w")
f.write('# LAMMPS data file \n\n')
f.write(str(cptatm)+' atoms\n')
f.write('\n')
f.write(str(int(2))+' atom types\n')
f.write('\n')
f.write(str(txlo)+' '+str(txhi)+' xlo xhi\n')
f.write(str(tylo)+' '+str(tyhi)+' ylo yhi\n')
f.write(str(tzlo)+' '+str(tzhi)+' zlo zhi\n')
f.write('\n')
f.write('Atoms\n')
f.write('\n')
for nlin in range(len(atoms)):
    newline = atoms[nlin]
    for col in range(len(newline)):
        if col < 3:
            f.write(str(int(newline[col]))+' ')
        else :
            f.write(str(newline[col])+' ')
    f.write('\n')
f.close()

