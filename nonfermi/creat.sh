#!/bin/bash
num=10
filename="nonfermi"
for (( i=1; i <= $num; i++))
do
    mkdir $filename$i
    cp skyrmionlattice $filename$i
    cp run.pbs $filename$i
    cp seed.in $filename$i
    #echo "0."$i > ./$filename$i/para.in
done
