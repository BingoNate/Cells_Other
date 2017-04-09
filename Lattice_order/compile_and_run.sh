#!/bin/bash

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Static_structure
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_1.0/fp_0.5/areak_10.0
point2=eps_1.0_fp_0.5_areak_10.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Glob_lattice_order
mkdir -p ${savebase}
savefile=${savebase}/Glob_lattice_order_${point2}.txt

cd ${path}
g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_global_lattice_order.cpp ${upath}/read_write.cpp -lhdf5 -o calc_glob_lat
./calc_glob_lat out.h5 ${savefile}
cd -

