## Water in a carbon nanotube using grand canonical Monte Carlo and molecular dynamics 

### Description

The first part of the run consists of a GCMC step. During this step, water molecules will be added inside the CNT. See "Understanding Molecular Simulation" by Daan Frenkel and Berend Smit for a description of GCMC algorithm. The second part is a MD step, during which the CNT is allowed to deform. 

![Algorithm schema](./WaterInCNT.jpeg)

### How to

Run the GenerateData.m using octave of matlab, or just use the data.lammps already generated. Then, run the input.lammps using LAMMPS.

### Output

The following video has been made with this code : https://www.youtube.com/watch?v=fIAmqMLPaZw

### See also

My LAMMPS tutorials website : https://lammpstutorials.github.io/


