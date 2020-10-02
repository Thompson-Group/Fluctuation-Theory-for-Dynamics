#!/usr/bin/env python
#grab_flucts.py is the program that parses each log file and grabs each of the energies/volumes/whatever thermodynamic data you want to weight the correlation fucntions. 
#   required input: either log.lammps or flucts.ener
#   produced output: *_init.out

import numpy as np
import sys
from read_input import user_input

input_file = str(sys.argv[1])
value = int(sys.argv[2])

inputparam = user_input('input_file')

filenames='file_names'

ivalname, icolnum = np.genfromtxt(input_file, usecols=(0,1), dtype=(str,str), unpack=True)
ival_list = np.array([])
if "elec" in ivalname[value]:
    sys.exit("Elec contribution calculated as the sum of the ewald and coul contributions later")

fileout=str(ivalname[value])+"_init.out"
outputfile=open(fileout,'w')
if inputparam.prog == 'LAMMPS':
    # Looks for log.lammps
    with open(filenames) as f:
        for l in f:
            filename="log.lammps"
            if float(l.rstrip()) < 501001000:
                filename='FILES/'+l.rstrip()+'/log.lammps'
            else:
                filename='FILES2/'+l.rstrip()+'/log.lammps'
            print("Grabbing File: %s" % filename)
            lookup='Step'
            lookup2='Loop time'
            with open(filename) as myFile:
                startskip=-1
                for num, line in enumerate(myFile,1):
                    if lookup in line:
                        startskip=num
                    if num == startskip+1:
                        ival=line.split()[int(icolnum[value])-1]
                        break
            outputfile.write("%s\n"%float(ival))
elif inputparam.prog == 'CP2K':
    # Looks for CP2K Output File
    with open(filenames) as f:
        for l in f:
            filename='FILES/'+l.rstrip()+'/free-1.ener'
            print("Grabbing File: %s" % filename)
            with open(filename) as myFile:
                myFile.readline()
                ival=myFile.readline().strip().split()[int(icolnum[value])-1]
                ival_list = np.append(ival_list, float(ival))
    fileout=str(ivalname[value])+"_init.out"
    np.savetxt(fileout, ival_list, fmt=['%.4f'])

            
else:
    print("%s is not a valid program" % inputparam.prog)

