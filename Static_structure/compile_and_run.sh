#!/bin/bash

set -e

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Static_structure
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_1.0/fp_0.5/areak_10.0/kappa_100.0
point2=eps_1.0_fp_0.5_areak_10.0_kapppa_100.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Static_struct
mkdir -p ${savebase}
savefile=${savebase}/Static_struct_${point2}.txt

cd ${path}
g++ -O3 -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/analyse_static_structure_log.cpp ${upath}/read_write.cpp -lhdf5 -lm -fopenmp -o static_struct
./static_struct out.h5 ${savefile}
cd -

