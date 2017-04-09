#!/bin/bash

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Pair_correlation
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_1.0/fp_0.5/areak_10.0
point2=eps_1.0_fp_0.5_areak_10.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Pair_corr
savefile=${savebase}/pair_corr_${point2}.txt

cd ${path}
g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calculate_pair_corr.cpp ${upath}/read_write.cpp -lhdf5 -o paircorr
./paircorr out.h5 ${savefile}
cd -

