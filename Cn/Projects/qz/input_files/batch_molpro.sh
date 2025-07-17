#!/bin/bash 
#mkdir ${PWD}/molpro_tmp

export TMPDIR="/DATA/molpro_tmp" 
export TMPDIR4="/DATA/molpro_tmp"

cd Molpro_CP 

for k in 356
do 
/usr/local/molpro/molprop_2015_1_linux_x86_64_i8/bin/molpro -n 32 $k.inp 
done 


