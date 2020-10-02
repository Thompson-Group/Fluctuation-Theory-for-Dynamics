#!/usr/bin/env python
#set_msd_calcs.py
""" 
This file creates the corr_calc.in file that is needed by the correlation function calculations to run them.

This file is used by msd_rot_calc among others.
"""

from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter
from read_input import user_input

# Reads the input file
inputparam = user_input('../../input_file')
fluctval = int(os.path.split(os.getcwd())[1])


if inputparam.timestep == 'FALSE':
    timestep=1.0
else:
    timestep = float(inputparam.timestep)
with open('mol.info') as f:
        mol_name = str(f.readline().strip())

if inputparam.prog == "LAMMPS":
    # Looks up the Step and Loop time keywords from log.lammps
    filename='log.lammps'
    lookup= 'Step'
    lookup2= 'Loop time'
    totlines=0
    startskip=0
    endskip=0
    with open(filename) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                startskip=num
            if lookup2 in line:
                endskip=num
            totlines=num
    endskip=totlines-endskip-4
    volu = np.genfromtxt(filename, skip_header=startskip,skip_footer=endskip,usecols=(11), unpack=True)
    vol = volu[0]
elif inputparam.prog == "CP2K":
    with open("in.nve.cp2k","r") as fi:
        id = 0
        for ln in fi:
            if ln.strip().startswith("ABC"):
                id=ln.strip().split()[-1]
        vol = float(id)**3
        print(vol)


# Writes the new msd calc input file with the right molecule name and volume
filepath='corr_calc.in'
f=open(filepath, 'w')

f.write("# Numeric File\n")
f.write("%s\n" % (fluctval))
f.write("# Number_of_Times Sep_of_Times:\n")
f.write("%s %s\n" % (inputparam.num_times, timestep))
f.write("# Volume\n")
f.write("%s\n" % (vol))
f.write("# Molecule Name\n")
f.write("%s\n" % (mol_name))
f.write("# Constraint (if applicable)\n")
f.write("%s\n" % (inputparam.constraint))
f.write("# Programs (1 for LAMMPS 2 for CP2K)\n")
if inputparam.prog == "LAMMPS":
    f.write("1\n")
elif inputparam.prog == "CP2K":
    f.write("2\n")
else:
    f.write("0\n")

f.close()
