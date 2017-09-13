#!/bin/bash
num=9
cd data
touch energy.out
cd ../
filename="nonfermi"
for (( i=0; i<=$num; i++))
do
    cat $filename$i/enenrgy.out >> data/energy.out
done
