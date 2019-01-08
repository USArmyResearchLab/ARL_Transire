#!/bin/bash

cd  NON_ENERGY_INPUTS/
for file in *.in;
do
    echo " "
    echo "============================"
    echo "Beginning test input : $file"
    echo "============================"
    echo " "
    transire.py -i "$file"
done
cd ..
