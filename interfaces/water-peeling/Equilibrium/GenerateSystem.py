#!/usr/bin/env python
# coding: utf-8

# In[17]:


from math import *

fichier = open("temp", "w")

# PARAMETERS 
typat = 4 #type d atomes
typbo = 1 #type de bond
typan = 1 #type d angle
nbbond = 0 #nombre bond
nbat=0 #nombre atome
nban=0 #nombre d angle

# Masses
mH = 1.00794
mO = 15.9994
mC = 12.0107

# Parametres wor water molecule geometry
alpha=104.5*pi/180
r=0.9584

# Parameter for graphene geometry
d=1.4 # airebo
beta=30*pi/180

# Charges
qO=-1.1128
qH=0.5564
qC=0

# Box definition
e=3.1589
sigma=3.374 #A
unit=2*sigma/sqrt(2)

c=6.696/2
L=6*unit

pos_wall=70
H2=L*5 # H*5

txlo=-34*d*cos(beta); txhi=-txlo
tylo=-3*(2*d*sin(beta)+2*d); tyhi=-tylo
tzlo=-30; tzhi=-tzlo
tzlo=-c/2-c-10
tzhi=int(pos_wall)+20

dX=6;
space=9;

# Creation systeme
fichier.write("# membrane \n \n")

# Box
fichier.write(str(txlo) + " " + str(txhi) + " xlo xhi \n")
fichier.write(str(tylo) + " " + str(tyhi) + " ylo yhi \n")
fichier.write(str(tzlo) + " " + str(tzhi) + " zlo zhi \n")
fichier.write("\n")

# Masses
fichier.write("Masses \n \n")
fichier.write("1 " + str(mC) + "\n")
fichier.write("2 " + str(mC) + "\n")
fichier.write("3 " + str(mO) + "\n")
fichier.write("4 " + str(mH) + "\n")
fichier.write("\n")

# Atoms
fichier.write("Atoms \n \n")

# Wall placement
nbat=0
nbcar=0
cptE=0

#### Graphene

for cpt, z in enumerate([-c*2, -c , 0, c]):
    x=txlo+(cpt%2)*d*cos(beta)/2
    y=tylo
    while y <= tyhi:
        while x <= txhi:
            x0=x
            y0=y
            z0=z
            if y0<tyhi and y0>=tylo and x0<txhi and x0>=txlo:
                fichier.write(str(nbcar+1) + " 0 " + "1 " + str(qC) 
                              + " " + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
            x0=x
            y0=y+d
            z0=z
            if y0<tyhi and y0>=tylo and x0<txhi and x0>=txlo:
                fichier.write(str(nbcar+1) + " 0 " + "1 " + str(qC) 
                              + " " + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n"); 
                nbat += 1
                nbcar += 1
            x0=x+d*cos(beta)
            y0=y+d*sin(beta)+d
            z0=z
            if y0<tyhi and y0>=tylo and x0<txhi and x0>=txlo:
                fichier.write(str(nbcar+1) + " 0 " + "1 " + str(qC) 
                              + " " + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
            x0=x+d*cos(beta)
            y0=y+d*sin(beta)+d*2
            z0=z
            if y0<tyhi and y0>=tylo and x0<txhi and x0>=txlo:
                fichier.write(str(nbcar+1) + " 0 " + "1 " + str(qC) 
                              + " " + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
            x=x+2*d*cos(beta)
        x=txlo+(cpt%2)*d*cos(beta)/2
        y=y+2*d*sin(beta)+2*d
        
#### Shearing walls
x=txlo+unit/4;
y=tylo+unit/4;
z=pos_wall+dX;
p=sqrt(2)*unit/2*0.7071067811865476;
nwall=0;
while y < tyhi:
    while x < txhi:
        while z <= pos_wall+unit*1.1+dX:
            x0=x
            y0=y
            z0=z+unit/2
            if x0<txhi and x0>=txlo and y0<tyhi and y0>=tylo:
                fichier.write(str(nbcar+1) + " 0 " + "2 " + " 0 "
                              + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
                nwall+=1
            x0=x;
            y0=y+unit/2;
            z0=z;
            if x0<txhi and x0>=txlo and y0<tyhi and y0>=tylo:
                fichier.write(str(nbcar+1) + " 0 " + "2 " + " 0 "
                              + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
                nwall+=1
            x0=x+unit/2
            y0=y
            z0=z
            if x0<txhi and x0>=txlo and y0<tyhi and y0>=tylo:
                fichier.write(str(nbcar+1) + " 0 " + "2 " + " 0 "
                              + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
                nwall+=1
            x0=x+unit/2
            y0=y+p
            z0=z+p
            if x0<txhi and x0>=txlo and y0<tyhi and y0>=tylo:
                fichier.write(str(nbcar+1) + " 0 " + "2 " + " 0 "
                              + str(x0) + " " + str(y0) + " " 
                              + str(z0) + " " + "\n")
                nbat += 1
                nbcar += 1
                nwall+=1
            z=z+unit;
        x=x+unit
        z=pos_wall+dX
    x=txlo+unit/4
    y=y+unit

print('Number of atom in the shearing wall : ' + str(nwall))

# placement des molecules d'eau
x=txlo+e/2
y=tylo+e/2
z=c/2+c+e

while x <= txhi-e/2:
    while y <= tyhi-e/2 :
        while z <= pos_wall-e/2 :
            if z**2 > space**2 :
                fichier.write(str(nbat+1) + " 1 " + "3 " + str(qO) 
                              + " " + str(x) + " " + str(y) + " " 
                              + str(z) + "\n")
                nbat += 1
                fichier.write(str(nbat+1) + " 1 " + "4 " + str(qH) 
                              + " " +  str(r+x)  + " " + str(y) + " " 
                              + str(z) + "\n") 
                nbat += 1
                fichier.write(str(nbat+1) + " 1 " + "4 " + str(qH) 
                              + " " + str(r*cos(alpha)+x) + " " + str(r*sin(alpha)+y) 
                              + " " + str(z) + "\n")
                nbat += 1
            z=z+e
        y=y+e
        z=c/2+c+e
    x=x+e
    y=tylo+e/2
    
fichier.write("\n");

cpth=(nbat-nbcar)/3
cptO=1

print('Number of water molecules : ' + str(cpth))

# Bonds
fichier.write("Bonds \n \n")
while cptE <= cpth-1 :		
    fichier.write(str(nbbond+1) + " 1 " + str(cptO+nbcar) + " " + str(cptO+1+nbcar) + "\n")
    fichier.write(str(nbbond+2) + " 1 " + str(cptO+nbcar) + " " + str(cptO+2+nbcar) + "\n")
    nbbond += 2
    cptO += 3
    cptE +=1

fichier.write("\n");

# Angles
fichier.write("Angles \n \n")
cptt=0;
while nban <= cpth-1 :
    fichier.write(str(nban+1) + " 1 " + str(cptt+2+nbcar) + " " 
                  + str(cptt+1+nbcar) + " " + str(cptt+3+nbcar) + "\n")
    nban +=1
    cptt += 3

fichier.write("\n")

# ecriture du fichier data 	
fichier.close()
source = open("temp", "r")
txt = source.read()
source.close()
import os
os.remove("temp")
fichier = open("data.lammps", "w")
# Introduction
fichier.write("#  Sampson"+ "\n"+ "\n")
fichier.write(str(nbat) + " " + "atoms"+ "\n")
fichier.write(str(nbbond) + " " + "bonds"+ "\n")
fichier.write(str(nban) + " " + "angles"+ "\n")
fichier.write("\n");
fichier.write(str(typat) + " " + "atom types"+ "\n")
fichier.write(str(typbo) + " " + "bond types"+ "\n")
fichier.write(str(typan) + " " + "angle types"+ "\n")
fichier.write("\n")
fichier.write(str(txt))
fichier.close()

