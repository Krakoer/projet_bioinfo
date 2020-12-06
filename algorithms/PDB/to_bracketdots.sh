#!/bin/bash


for f in pdb_files/*.pdb; 
do 
    filename="${f:10:4}"
    if test -f "outputs/$filename.dbn"; then
        echo "File $f already processed";
    else
        echo "Processing $f file.."; 
        ./x3dna-dssr "-i=$f";
        cp dssr-2ndstrs.dbn "outputs/$filename.dbn";
        rm dssr*;
    fi
done

# f="pdb_files/5LMQ.pdb";
# filename="${f:10:4}"
# echo "Processing $f file.."; 
# ./x3dna-dssr "-i=$f";
# cp dssr-2ndstrs.dbn "outputs/$filename.dbn";
# rm dssr*;