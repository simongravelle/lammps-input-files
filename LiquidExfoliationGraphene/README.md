## Liquid exfoliation of a multilayer nanographene: out-of-equilibrium molecular dynamic simulations

### Description

The simulation consists of a nanoparticle made of multiple graphene layers in a box of water. A linear shear is imposed thanks to two moving walls.  With a shear rate of 1e10 1/s or higher, the layers separate from each other and the particle is exfoliated into multiple single layer particles.

![Algorithm schema](./LiquideExfoliationGraphene.jpeg)

### How to

From the 'Construction' folder, run the GenerateData.m script using Matlab (Octave should also work with a few modification). From the 'Equilibration' folder, run the input.lammps script using LAMMPS. Finally, run the input.lammps script within the 'Run' folder using LAMMPS. 

### Output

The following video has been made with this code : https://www.youtube.com/watch?v=GALFLXkUEAU

### See also

My LAMMPS tutorials website : https://lammpstutorials.github.io/

Our article on the subject : https://doi.org/10.1063/1.5141515


