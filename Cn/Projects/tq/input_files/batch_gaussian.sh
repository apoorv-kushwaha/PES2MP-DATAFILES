#!/bin/bash 

# Load Gaussian with scratch dir information
. /home/user/Desktop/g16_lin/gaussian_vars_16.sh

cd Gaussian_CP

for k in {0..1424}
do 
g16 $k.gjf
done 


