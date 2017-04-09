#!/bin/sh

for P in 0.3 0.4 0.5 0.6 0.7; do
for E in 10 12 14; do
for F in 1 2; do

path=phi${P}eps${E}fm${F}
echo ${path}
cp plot_neigh.py ../Graphs/${path}/
cd ../Graphs/${path}
anaconda plot_neigh.py 50000 600 ./ num_neigh_dist
anaconda plot_neigh.py 50000 600 ./ num_neigh_per_frame
anaconda plot_neigh.py 50000 600 ./ cluster_size_dist
cd ../../Analysis

done
done 
done

