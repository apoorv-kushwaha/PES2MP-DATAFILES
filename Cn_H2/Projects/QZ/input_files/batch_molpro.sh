#!/bin/bash 
mkdir ${PWD}/molpro_tmp

export TMPDIR="/${PWD}/molpro_tmp" 
export TMPDIR4="/${PWD}/molpro_tmp"

cd Molpro_CP 

for k in {0..8399}
do 
/usr/local/molpro/molprop_2015_1_darwin_x86_64_i8/bin/molpro -n 2 $k.inp 
done 


