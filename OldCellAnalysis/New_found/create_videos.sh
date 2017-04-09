#!/bin/sh

for P in 0.3; do
for E in 10 14; do
for F in 1 2; do

path=phi${P}eps${E}fm${F}

echo Analysis for ${path}

cp ../Folders/${path}/variables.cpp .

echo Field calculation
icpc variables.cpp field.cpp -o postofi -O3 -fp-model precise -no-multibyte-chars
./postofi ../Data/${path}/Data1/coords.xyz

cp img/dat/fastest_cells_1.txt img/dat/fastest_cells_0.txt
cp img/dat/cells_1.txt img/dat/cells_0.txt

echo Pressure calculation
icpc variables.cpp pressure.cpp -o postopr -O3 -fp-model precise -no-multibyte-chars
./postopr ../Data/${path}/Data1/coords.xyz

cp img/dat/local_vir_stress_1.txt img/dat/local_vir_stress_0.txt
cp img/dat/local_pressure_1.txt img/dat/local_pressure_0.txt

anaconda video_analysis.py 50000 600 ../Data/${path}/Data1/cell.xyz
cd img
anaconda st2_video.py
mv video-3.mpeg ${path}.mpeg
mv *.mpeg ../../Graphs/${path}/.
cd ..


done
done
done

