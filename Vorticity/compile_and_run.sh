#!/bin/bash

sname=Vorticity
dname=Vorticity_field
dname2=Enstrophy
fname=calc_vorticity

spath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/${sname}
upath=/usr/users/iff_th2/duman/Cells_in_LAMMPS/Scripts/Utility
point1=eps_1.0/fp_5.0/areak_10.0
point2=eps_1.0_fp_5.0_areak_10.0

dpath=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8
path=${dpath}/${point1}

savebase=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/${dname}
mkdir -p ${savebase}
savefile=${savebase}/${dname}_${point2}.txt

savebase2=/usr/users/iff_th2/duman/Cells_in_LAMMPS/DATA/${dname2}
mkdir -p ${savebase2}
savefile2=${savebase2}/${dname2}_${point2}.txt

cd ${path}
/usr/local/gcc/bin/g++ -std=c++11 -O3 -lm -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/${fname}.cpp ${upath}/read_write.cpp -lhdf5 -o ${fname}
./${fname} out.h5 ${savefile} ${savefile2}
cd -

