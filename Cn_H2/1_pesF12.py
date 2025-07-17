#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
#-------------------------    Important Options     ---------------------------#
################################################################################
#------------------------------------------------------------------------------#
# Enter parameter for RR1 : Enter experimental/theoretical bond length(s)------#
RR1_atoms        = ['C','C'] # Enter rigid rotor (RR1) atoms from one end
RR1_bond_len     = [1.2425]  # Enter bond lengths from the same end (see manual)
# Enter parameter for second (collider) atom ----------------------------------#
RR2_atoms        = ['H','H'] # Enter atom(s) for RR2
RR2_bond_len     = [0.7667]  # Enter bond lengths from the same end for RR2
#------------------------------------------------------------------------------#
# Enter charge & multiplicity (Do NOT change order) ---------------------------#
Charge           = [0, 0, 0] # charge of [atom 1, atom 2, whole system]
Multiplicity     = [1, 1, 1] # multiplicity of [atom 1, atom 2, whole system]
#------------------------------------------------------------------------------#

################################################################################
#-----------------------    2D/4D Plot Options     ----------------------------#
################################################################################
# 4D angular [phi, theta2, theta1]  coordinates (ab initio calculation) -------#
phi    = [0,  91, 30]    # phi    0 to 90 with step size of 30
theta2 = [0,  91, 30]    # theta2 0 to 90 with step size of 30
theta1 = [0, 181, 30]    # theta1 0 to 90 with step size of 30
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

################################################################################
################################################################################
#----------      Codes that probably don't need changes         ---------------#
################################################################################
################################################################################
#------------------------------------------------------------------------------#
# some necessary pieces of code for running job and grabbing project name -----#
Create_PES_input = True       # Set JOB Type(s) (True/False)
#------------------------------------------------------------------------------#
import os        # Gettig project name from GUI interface ---------------------#
Proj_name        =  os.getenv("Proj_name", "default_project_name") # Auto-set
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#---------------------- Advanced PES/Plot Options  ----------------------------#
#------------------------------------------------------------------------------#
# Radial coordinates (enter as many range and step sizes as desired) ----------#
R_i   = [2.5,  6.0,  9.0, 15.0,  50.0] # initial R in Ang. (eg. 0.1 Ang.)
R_f   = [6.0,  9.0, 15.0, 25.1, 100.1] # final R - R_stp (eg 1.95 Ang.)
R_stp = [0.1,  0.2,  0.5,  1.0,  50.0] # step size for R_i-R_f
# Loops run from initial to penultimate value (eg. loop 2-6 goes till 5.9 Ang.)#
#------------------------------------------------------------------------------#
#--------------------  Use Isotope (Optional)  --------------------------------#
#Use_isotope = True # only for calculating isotopic dependent COM (PES is same)!
#RR1_isotope_mass = [0, 13.00335483507]     # Enter isotope mass  (0 = default)
#RR2_isotope_mass = [0,  3.0160492779 ]     # Enter isotope mass  (0 = default)
#------------------------------------------------------------------------------#

Create_MOLPRO_CP_input_files      = True       # MOLPRO(CP)

mem_m      = '1 g'                    # memory per thread/processor (in words)
run_method = '''{{rhf;rccsd(t)-f12b,SCALE_TRIP=1}}'''   # molpro function to declare method
basis_cp   = 'AVDZ'                   # basis set for cp calculation
################################################################################
################################################################################
