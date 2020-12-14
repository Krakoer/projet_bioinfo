#!/bin/bash

for f in ../../RNA_files/pdb_files/*.pdb; 
do 
    filename="${f:26:4}"

    if test -f "../../RNA_files/json/$filename.json"; then
        echo "File $f already processed";
    else
        echo "Processing $f file.."; 
        cmd=" -i=$f -o=../../RNA_files/json/$filename.json --json";
        ./x3dna-dssr $cmd;
        rm dssr*;
    fi
done