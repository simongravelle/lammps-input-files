
# coding: utf-8

import numpy as np
import random
import copy

Na = 6.022e23 #constants.Avogadro
Mh2o = 0.018053 # kg/mol - water

# choose the initial box dimention
dnacl = 2.84
nx = 12
ny = 12
nz = 8
dw = 4.5
layer0 = nz*dnacl
layer = layer0
h = 40
Lx = nx*dnacl
Ly = ny*dnacl
Lz = layer0 + h
# initial NaCl mesh
Nnacl = (Lx/dnacl)*(Ly/dnacl)

multi = 0.3
N = np.int64(Nnacl*multi)

cptH2O = 0
nCl = 0
nNa = 0

Natomtypes = 4
Nbondtypes = 1
Nangletypes = 1

txlo, txhi = -Lx/2, Lx/2
tylo, tyhi = -Ly/2, Ly/2
tzlo, tzhi = -Lz/2, Lz/2
attemps = 0


# Load NaCl positions for the wall crystal structure and water
wallNaCl = np.zeros((10000,7))
file1 = open('NaCl/position.txt', 'r')
Lines = file1.readlines()
count = 0
for line in Lines:
    wallNaCl[count]=line.strip().split(' ')
    count += 1
wallNaCl = wallNaCl[0:count]

Ph2o = np.loadtxt('../ff/H2O_TIP4P2005/position.dat')
Bh2o = np.loadtxt('../ff/H2O_TIP4P2005/bond.dat')
Ah2o = np.loadtxt('../ff/H2O_TIP4P2005/angle.dat')

while cptH2O+nNa+nCl < N: # if the initial distance between the wall is too small to acomodate the molecules, increase h
    if attemps>0:
        print('Increasing the distance between the walls')
        h += dw
        print('new h = '+str(h)+ ' Angstroms')
        Lz = layer0 + h 
        tzlo = -Lz/2 
        tzhi = Lz/2
        
    cptatom = 0
    cptbond = 0
    cptangle = 0
    cptmol = 0
    cptNa = 0
    cptH2O = 0
    cptCl = 0
    nCl = 0
    nNa = 0
    
    atoms = np.zeros((10000,7))
    bonds = np.zeros((10000,4))
    angles = np.zeros((10000,5))
    
    # replicate the initial structure
    wallNaClrep = []
    wallNaClrep = copy.deepcopy(wallNaCl)
    for xx in np.arange(txlo+dnacl/2,txhi,2*dnacl):
        for yy in np.arange(tylo+dnacl/2,tyhi,2*dnacl):
            for zz in np.arange(tzlo+dnacl/2,tzlo+layer0,2*dnacl):
                wallNaClrep = np.append(wallNaClrep,wallNaCl+[0,0,0,0,xx,yy,zz], axis=0)
    wallNaClrep = wallNaClrep[8:]
    assert len(wallNaClrep.T[6])>0
    wallNaClrep = wallNaClrep[wallNaClrep.T[6] < np.max(wallNaClrep.T[6])]
    assert len(wallNaClrep[wallNaClrep.T[2]==1]) == len(wallNaClrep[wallNaClrep.T[2]==2])

    for n in range(len(wallNaClrep)):
        wallNaClrep[n,0] = np.int64(n+1)
        
    minwall = np.min(wallNaClrep.T[6])
    maxwall = np.max(wallNaClrep.T[6])
    layer = maxwall-minwall
    wallNaClrep.T[6] -= (minwall-tzlo)
    wallNaClrep.T[6] -= layer/2
    
    # build the wall
    cptmol += 1
    for m in wallNaClrep:
        atoms[cptatom] = cptatom+1, cptmol, m[2]+2, m[3], m[4], m[5], m[6]
        cptatom += 1
    cptwall = cptatom
  
    # replicate all the atoms for neighbor search with periodic boundary condition
    # note : this could be improved using the search neighbor function of MDAnalysis  
    atomsrep = copy.deepcopy(atoms[0:cptwall])
    for xx in [-Lx,0,Lx]:
        for yy in [-Ly,0,Ly]:
            for zz in [-Lz,0,Lz]:
                atomsrep = np.append(atomsrep,atoms[0:cptwall]+[0,0,0,0,xx,yy,zz], axis=0)  
    atomsrep = atomsrep[cptwall:]
    
    for z in np.arange(tzlo+dw/2,tzhi-dw/4,dw):
        for x in np.arange(txlo+dw/2,txhi-dw/4,dw):
            for y in np.arange(tylo+dw/2,tyhi-dw/4,dw):
            
                d = np.sqrt((x-atomsrep.T[4])**2+(y-atomsrep.T[5])**2+(z-atomsrep.T[6])**2)
                rdm = random.random()
                
                if (cptH2O < N) & (np.min(d) > dw) & (rdm<0.2):
                    
                    cptmol += 1
                    cptH2O += 1
                    for m in Bh2o:
                        bonds[cptbond] = cptbond+1, m[1], m[2]+cptatom, m[3]+cptatom
                        cptbond += 1

                    m = Ah2o
                    angles[cptangle] = cptangle+1, m[1], m[2]+cptatom, m[3]+cptatom, m[4]+cptatom
                    cptangle += 1 

                    for m in Ph2o:
                        atoms[cptatom] = cptatom+1, cptmol, m[2], m[3], m[4]+x, m[5]+y, m[6]+z
                        cptatom += 1 
    attemps += 1
atoms = atoms[0:cptatom]       
bonds = bonds[0:cptbond]    
angles = angles[0:cptangle]

# write LAMMPS data file
f = open("data.lammps", "w")
f.write('# LAMMPS data file \n\n')
f.write(str(cptatom)+' atoms\n')
f.write(str(cptbond)+' bonds\n')
f.write(str(cptangle)+' angles\n')
f.write('\n')
f.write(str(Natomtypes)+' atom types\n')
f.write(str(Nbondtypes)+' bond types\n')
f.write(str(Nangletypes)+' angle types\n')
f.write('\n')
f.write(str(-Lx/2)+' '+str(Lx/2)+' xlo xhi\n')
f.write(str(-Ly/2)+' '+str(Ly/2)+' ylo yhi\n')
f.write(str(-Lz/2)+' '+str(Lz/2)+' zlo zhi\n')
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
f.write('\n')
f.write('Bonds\n')
f.write('\n')
for nlin in range(len(bonds)):
    newline = bonds[nlin]
    for col in range(len(newline)):
        f.write(str(int(newline[col]))+' ')
    f.write('\n')
f.write('\n')
f.write('Angles\n')
f.write('\n')
for nlin in range(len(angles)):
    newline = angles[nlin]
    for col in range(len(newline)):
        f.write(str(int(newline[col]))+' ')
    f.write('\n')
f.close()

