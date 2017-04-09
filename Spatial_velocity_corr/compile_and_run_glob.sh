#!/bin/bash

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Spatial_velocity_corr
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_1.0/fp_0.5/areak_10.0
point2=eps_1.0_fp_0.5_areak_10.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Sp_velocity_corr
savebase2=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Vel_corr_length
mkdir -p ${savebase}
mkdir -p ${savebase2}
savefile=${savebase}/Sp_velocity_corr_${point2}
savefile2=${savebase2}/Vel_corr_length_${point2}.txt

cd ${path}
/usr/local/gcc/bin/g++ -std=c++11 -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_vel_corr_length.cpp ${upath}/read_write.cpp -lhdf5 -o calc_vel_corr_length
./calc_vel_corr_length out.h5 ${savefile} ${savefile2}
cd -

