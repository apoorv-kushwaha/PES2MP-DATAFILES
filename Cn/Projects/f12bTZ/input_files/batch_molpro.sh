#!/bin/bash 
mkdir ${PWD}/molpro_tmp

export TMPDIR="/${PWD}/molpro_tmp" 
export TMPDIR4="/${PWD}/molpro_tmp"

cd Molpro_CP 

for k in {0..1424}
do 
/usr/local/molpro/molprop_2015_1_darwin_x86_64_i8/bin/molpro -n 10 $k.inp 
done 


