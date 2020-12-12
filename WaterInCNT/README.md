Grand Canonical Monte Carlo simulation + Molecular dynamics simulation of water molecule in an effectively infinite carbon nanotube.

Step 1: Run the GenerateData.m using octave of matlab, or just use the data.lammps already generated.

Step 2: Run the input.lammps using LAMMPS.

The first part of the run consists of a GCMC step. During this step, water molecules will be added inside the CNT. See "Understanding Molecular Simulation" by Daan Frenkel and Berend Smit for a description of GCMC algorithm. 

The second part is a MD step, during which the CNT is allowed to deform. 
