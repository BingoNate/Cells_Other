#!/bin/sh

cd ..
mkdir Graphs

for P in 0.3; do
for E in 6 10 22; do
for F in 1; do

path=phi${P}eps${E}fm${F}
echo ${path}
mkdir Graphs/${path}
cp Folders/${path}/variables.cpp Analysis/variables2.cpp
cd Analysis
icpc variables2.cpp bondlife.cpp -o pt -O3 -fp-model precise -no-multibyte-chars
./pt ../Data/${path}/Data1/coords.xyz
#icpc variables2.cpp main.cpp -o ma -O3 -fp-model precise -no-multibyte-chars
#./ma ../Data/${path}/Data1/coords.xyz
#icpc variables2.cpp bond_corr.cpp -o pt -O3 -fp-model precise -no-multibyte-chars
#./pt ../Data/${path}/Data1/coords.xyz
#icpc variables2.cpp velocity_corr.cpp -o pt -O3 -fp-model precise -no-multibyte-chars
#./pt ../Data/${path}/Data1/coords.xyz
#icpc variables2.cpp density_corr.cpp -o pt -O3 -fp-model precise -no-multibyte-chars
#./pt ../Data/${path}/Data1/coords.xyz
#icpc variables2.cpp msd.cpp -o pt -O3 -fp-model precise -no-multibyte-chars
#./pt ../Data/${path}/Data1/coords.xyz
#icpc variables2.cpp virial.cpp -o pt -O3 -fp-model precise -no-multibyte-chars
#./pt ../Data/${path}/Data1/coords.xyz

#anaconda plot_postData.py
ls | grep txt | grep -v avg | grep -v max | grep -v susc | grep -v glob | xargs -I blah mv blah ../Graphs/${path}/.
ls | grep png | xargs -I blah mv blah ../Graphs/${path}/.
cd ..

done
done
done

