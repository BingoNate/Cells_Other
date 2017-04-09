#!/bin/bash

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Spatial_velocity_corr
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_20.0/fp_1.0/areak_1.0
point2=eps_20.0_fp_1.0_areak_1.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/Sp_velocity_corr
mkdir -p ${savebase}
savefile=${savebase}/Sp_velocity_corr_${point2}.txt

cd ${path}
g++ -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calc_sp_velocity_corr.cpp ${upath}/read_write.cpp -lhdf5 -o calc_sp_vel_corr
./calc_sp_vel_corr out.h5 ${savefile}
cd -

