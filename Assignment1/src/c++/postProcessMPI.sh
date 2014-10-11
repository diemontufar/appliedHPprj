#!/bin/bash
rm -f Assignment1_Files_MPI/*
for i in $(ls Shallow_Water_MPI.csv.*); do
    awk -F, '{print $2","$3","$4 >> "Assignment1_Files_MPI/Shallow_Water_MPI."$1".csv.tmp"}' $i
done
cd Assignment1_Files_MPI
for i in $(ls Shallow_Water_MPI*.csv.tmp); do
    name=$(echo $i | cut -d "." -f 1,2)
    sort -t\, -n -k 1,1n -k 2,2n -k 3,3n -k 4,4n $i > $name".csv"
done
cd ..
rm -f Assignment1_Files_MPI/*.csv.tmp
