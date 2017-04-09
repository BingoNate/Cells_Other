#!/bin/bash

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Intermediate_scattering
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_1.0/fp_0.5/areak_10.0
point2=eps_1.0_fp_0.5_areak_10.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Inter_scattering
mkdir -p ${savebase}
savefile=${savebase}/Inter_scattering_${point2}.txt

cd ${path}
g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calculate_inter_scattering.cpp ${upath}/read_write.cpp -lhdf5 -fopenmp -o calc_inter_scattering
./calc_inter_scattering out.h5 ${savefile}
cd -

