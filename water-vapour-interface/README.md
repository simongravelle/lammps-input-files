## Measure the surface tension of a water slab

### Description

Molecular dynamics simulation of a slab of water in contact with vacuum. The water model is TIP4P/2005. The surface tension is measured from the difference between the stresses normal and parallel to the liquid-vapour interface. 

![Algorithm schema](./water-vapour.png)

### How to

Run the make-data-file.m using Octave or Matlab to create the lammps data file. Then, execute the input.lammps file using LAMMPS. Visualise the dump file using  VMD, and extract the temperature and surface tension from the two data files. If you are new to LAMMPS and VMD, you can find [tutorials and instructions here](https://lammpstutorials.github.io/).

### Output

This [video](https://www.youtube.com/watch?v=l_APjA5_wZc) has been made with this script.

### Contact

Feel free to contact me by email if you have inquiries. You can find contact details on my [personal page](https://simongravelle.github.io/).
