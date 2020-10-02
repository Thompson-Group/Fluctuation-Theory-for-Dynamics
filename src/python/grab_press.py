#!/usr/bin/env python
"""
This is a code that pulls the pressures from the log file
"""
import numpy as np

it=0
with open ("log.lammps","r") as log:
    for line in log:
        it+=1
        if 'Pxx Pyy Pzz' in line:
            start=it
            
        if 'Loop time' in line:
            end=it
log.close()
P = np.genfromtxt('log.lammps', skip_header=start, max_rows=int(end-start-1), usecols=(13,14,15,16,17,18))

np.savetxt('pressures_out.log', P, fmt='%f8')


