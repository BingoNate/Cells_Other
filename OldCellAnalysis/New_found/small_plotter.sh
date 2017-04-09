#!/bin/sh

for P in 0.3 0.4 0.5 0.6 0.7; do
for E in 10 12 14; do
for F in 1 2; do

path=phi${P}eps${E}fm${F}
echo ${path}
cp plot_postData.py ../Graphs/${path}/
cd ../Graphs/${path}
anaconda plot_postData.py 50000 600 cell.xyz 
cd ../../Analysis

done
done 
done

