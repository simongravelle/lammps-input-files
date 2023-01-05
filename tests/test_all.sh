#!/bin/bash

set -e

lmp=/work/sgravelle/Softwares/lammps-3Nov2022/src/lmp_mpi
test=input_test.lammps

cd ..

cd solid-NaCl-under-compression/
    jupyter-nbconvert --to script make-LAMMPS-data.ipynb
    python3 make-LAMMPS-data.py
    cp input.lammps ${test}
    sed -i ' {n;/run / {s/run/run 2 #/}}' ${test}
    mpirun -np 4 ${lmp} -in ${test}
    rm dump.lammpstrj ${test} log.lammps make-LAMMPS-data.py NaCl.data
cd ..

cd liquid-mixture-PEG-water-ethanol/
    jupyter-nbconvert --to script create_mixture.ipynb
    python3 create_mixture.py
    cp input.lammps ${test}
    sed -i  ' {n;/run / {s/run/run 100 #/}}' ${test}
    mpirun -np 4 ${lmp} -in ${test}
    rm log.lammps dump.* mixture.gro mixture.data ${test} create_mixture.py peg.data
cd ..