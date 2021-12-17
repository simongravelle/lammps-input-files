#!/bin/sh 
#SBATCH --job-name=BulkTip42005nve
#SBATCH --mail-user=sgravelle@icp.uni-stuttgart.de
#SBATCH --mail-type=FAIL
#SBATCH -n 8
#SBATCH --time=48:00:0
#SBATCH --mem=1000

module load gcc
module load openmpi

srun -n 8 /home/sgravelle/mylammps/src/lmp_mpi -in npt.lammps
srun -n 8 /home/sgravelle/mylammps/src/lmp_mpi -in run.lammps
