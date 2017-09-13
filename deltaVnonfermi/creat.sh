#!/bin/bash
num=9
filename="nonfermi"
for (( i=0; i <= $num; i++))
do
    mkdir $filename$i
    cp skyrmionlattice $filename$i
    cp run.pbs $filename$i
    cp seed.in $filename$i
    echo "0.00"$i > ./$filename$i/B_z.in
done
