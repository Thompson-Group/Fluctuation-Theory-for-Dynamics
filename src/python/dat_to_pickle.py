#!/usr/bin/env python
import numpy as np
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-corr', default="msd", type=str, help='Filename without the ".dat"')
args = parser.parse_args()

corr = args.corr

with open('corr_calc.in','r') as f:
    lines=f.readlines()
    numeric = str(int(lines[1].strip()))
    mol = str(lines[7].strip())

base = corr + "_" + numeric + "_" + mol
data = np.genfromtxt(base + ".dat", usecols=1)
pickle.dump(data, open(base+".pckl",'wb'))
