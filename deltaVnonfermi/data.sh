#!/bin/bash
num=9
filename="nonfermi"
for (( i=0; i<=$num; i++))
do
    cp $filename$i/lattice.out data/lattice_$i.out
done
