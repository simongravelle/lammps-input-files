import numpy as np
import random
import copy

Na = 6.022e23 #constants.Avogadro
Mh2o = 0.018053 # kg/mol - water

dw = 3.2 		# initial inter-molecule distances (A)
Ntot = 2000		# desired number of species (ion + water)

c = 5
nion = c*Ntot*Mh2o/(2*(1+Mh2o*c)) # desired number for each ion
nwater = Ntot - 2*nion

# Types of atoms, bonds, and angles
Tatom = 4
Tbond = 1
Tangle = 1

# System dimensions
Lx, Ly, Lz = 25, 25, 25 # initial box size (A)
txlo, txhi = -Lx/2, Lx/2
tylo, tyhi = -Ly/2, Ly/2
tzlo, tzhi = -Lz/2, Lz/2

# Load water molecule
TypH2O = ['OW', 'HW1', 'HW2']
Ph2o = np.loadtxt('../../../shared/H2O_TIP4P2005/Position.dat')
Bh2o = np.loadtxt('../../../shared/H2O_TIP4P2005/Bond.dat')
Ah2o = np.loadtxt('../../../shared/H2O_TIP4P2005/Angle.dat')

cptH2O = 0
cptNa = 0
cptCl = 0
attemps = 0
while cptH2O+cptNa+cptCl < Nh2O:
    if attemps>0:
        Lx += dw/2
        Ly += dw/2
        Lz += dw/2
        txlo, txhi = -Lx/2, Lx/2
        tylo, tyhi = -Ly/2, Ly/2
        tzlo, tzhi = -Lz/2, Lz/2
        print('Increasing the box size, new Lx = '+str(Lx)+' A')
        
    cptatom = 0
    cptbond = 0
    cptangle = 0
    cptmol = 0
    cptH2O = 0
    cptCl = 0
    cptNa = 0
    cptCl = 0
    
    # allocate memory
    atoms = np.zeros((100000,7))
    bonds = np.zeros((100000,4))
    angles = np.zeros((100000,5))
    
    
    
    
    for z in np.arange(tzlo+dw/2,tzhi-dw/2,dw):
        for x in np.arange(txlo+dw/2,txhi-dw/2,dw):
            for y in np.arange(tylo+dw/2,tyhi-dw/2,dw):
                if cptH2O < Nh2O:
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

# write data.lammps
f = open("data.lammps", "w")
f.write('# LAMMPS data file \n\n')
f.write(str(cptatom)+' atoms\n')
f.write(str(cptbond)+' bonds\n')
f.write(str(cptangle)+' angles\n')
f.write('\n')
f.write(str(Tatom)+' atom types\n')
f.write(str(Tbond)+' bond types\n')
f.write(str(Tangle)+' angle types\n')
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
