#!/bin/bash
num=10
filename=nonfermi
for (( i=1; i<= $num; i++ ))
do
    cp $filename$i/lattice.out data/lattice_$i.out
done
