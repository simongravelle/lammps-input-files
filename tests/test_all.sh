#!/bin/bash

set -e

lmp=/work/sgravelle/Softwares/lammps-3Nov2022/src/lmp_mpi
test=input_test.lammps

cd ../liquid-mixture-PEG-water-ethanol/
    jupyter-nbconvert --to script create_mixture.ipynb
    python3 create_mixture.py
    cp input.lammps input_test.lammps
    sed -i  ' {n;/run / {s/run/run 100 #/}}' ${test}
    mpirun -np 4 ${lmp} -in ${test}
    rm log.lammps dump.* mixture.gro input_test.lammps create_mixture.py peg.data
