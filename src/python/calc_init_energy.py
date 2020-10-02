#!/usr/bin/env python
"""
This is a python code to calculate specific interactions and tabulate their interactions
    General Outline:
        1) read input file for pair coeffs
        2) read data file for types and molecules
        3) calculate interaction matrix
        4) caculate lennard jones potential
        5) calculate short range coulomb
"""

def read_paircoeffs(file):
    tmptype = []
    tmpeps =[]
    tmpsig = []
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            if "pair_coeff" in line:
                tmptype.append(int(line.split()[1]))
                tmpeps.append(float(line.split()[3]))
                tmpsig.append(float(line.split()[4]))
    eps = []
    sig = []
    for ty in tmptype:
        eps.append(0)
        sig.append(0)

    for i in range(len(tmptype)):
        eps[tmptype[i]-1] = tmpeps[i]
        sig[tmptype[i]-1] = tmpsig[i]
    return eps, sig

def pbc_dist(r1,r2,L):
    drij = r1 - r2
    drij = drij - L*np.round(drij/L)
    dr = np.sqrt(drij[0]**2. + drij[1]**2. + drij[2]**2.)
    return dr

def lb_mix(e1,e2,s1,s2):
    e = np.sqrt(e1*e2)
    s = (s1+s2)/2
    return e,s

def calc_LJ(dr, epsilon, sigma):
    dsr6 = (sigma/dr)**6.
    dsr12 = dsr6**2
    return 4*epsilon*(dsr12-dsr6)

def calc_COUL(dr, q1, q2):
    C = 332.0739
    return C*q1*q2/dr

def read_frame():
    return

def calc_totE():
    return

def write_E():
    return



if __name__ == "__main__":
    import sys
    import numpy as np
    
    inpfile = str(sys.argv[1])
    datafile = str(sys.argv[2])
    
    read_paircoeffs(inpfile)
    r1 = np.asarray([4,2,1])
    r2 = np.asarray([24,2,0])
    L = np.asarray([100.,100.0,100.0])
    calc_LJ(r1,r2,1,1,0.2,2.4,L)
    
