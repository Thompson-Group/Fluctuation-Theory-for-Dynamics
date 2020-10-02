#!/usr/bin/env python
import numpy as np
import argparse
from scipy import stats
from read_input import user_input

parser = argparse.ArgumentParser()
parser.add_argument('-etype', default='flucts.special.inp',type=str, help='What is the energy type that you are interested in?')
parser.add_argument('-corr_name', default='crp',type=str, help='What is the base corr_name that you are interested in?')
parser.add_argument('-dcorr_name', default='dcrp',type=str, help='What is the base dcorr_name that you are interested in?')
parser.add_argument('-mol_name', default='water',type=str, help='What is the molecule name?')

args = parser.parse_args()

fname_etype=args.etype
corr_name=args.corr_name
dcorr_name=args.dcorr_name
mol_name=args.mol_name


# read the input file
inputparam = user_input("input_file")

t_val = stats.t.ppf(0.975,inputparam.nblocks-1)/np.sqrt(inputparam.nblocks)

# read in the time

time = np.genfromtxt('real_time.dat', usecols=0)
time = [0] + time
time = time[:inputparam.num_times]

# read in other energies
coulen=np.genfromtxt('coul_init.out',usecols=0,unpack=True)
kspen=np.genfromtxt('ewald_init.out',usecols=0,unpack=True)

ljen=np.genfromtxt('lj_init.out', usecols=0, unpack=True)
elecen = np.add(coulen,kspen)

# read in file names
fnames = np.genfromtxt('file_names', dtype=str, unpack=True)

fcab=["FILES/"+str(s)+"/"+corr_name+"_"+str(s)+"_"+mol_name+".dat" for s in fnames]

cab={}
cab["cab"], cab["blcab"] = [],[]
tcab=[]
num_configs = int((inputparam.end_config+1000-inputparam.start_config)/inputparam.sep_config)
for i in range(num_configs):
    tcab.append(np.genfromtxt(fcab[i],usecols=(1),unpack=True))
cab["cab"]=np.average(tcab, axis=0)
    
for block in range(inputparam.nblocks):
    cab["blcab"]=np.average(np.split(np.array(tcab),inputparam.nblocks)[block],axis=0)


# read in energy
etypes=np.genfromtxt(fname_etype, usecols=(0),dtype=str, unpack=True)
energy={"reducedlj":ljen, "reducedelec":elecen}
aven={}
outcab={}
for etype in etypes:
    print("Now working on %s" % etype)
    # Specify corr file names
    print(dcorr_name, mol_name, etype)
    dfcab=["FILES/"+str(s)+"/"+etype+dcorr_name+"_"+str(s)+"_"+mol_name+".dat" for s in fnames]
    # Read energy and calculate averages
    energy[etype]=np.genfromtxt(etype+'_init.out', usecols=(0), unpack=True)
    if "LJ" in etype:
        energy["reducedlj"]=np.subtract(energy["reducedlj"],energy[etype])
    if "C" in etype:
        energy["reducedelec"]=np.subtract(energy["reducedelec"],energy[etype])

    aven[etype]=np.average(energy[etype])
    aven["bl"+etype],cab["blav"+etype]=[],[]
    # Read Corr Files
    cab[etype]=[]
    for i in range(num_configs):
        cab[etype].append(np.genfromtxt(dfcab[i],usecols=1, unpack=True))
    cab["av"+etype]=np.average(cab[etype],axis=0)
    for block in range(inputparam.nblocks):
        aven["bl"+etype].append(np.average(np.split(energy[etype],inputparam.nblocks)[block]))
        cab["blav"+etype].append(np.average(np.split(np.array(cab[etype]),inputparam.nblocks)[block],axis=0))

    # Calculate Derivative
    outcab["d"+etype]=cab["av"+etype]-np.multiply(aven[etype],cab["cab"])
    outcab["bld"+etype]=[]
    for block in range(inputparam.nblocks):
        outcab["bld"+etype].append(cab["blav"+etype][block]-np.multiply(aven["bl"+etype][block],cab["blcab"][block]))
    outcab["errd"+etype] = np.std(outcab["bld"+etype],axis=0)*t_val

np.savetxt("reducedlj_init.out", np.c_[energy["reducedlj"]])
np.savetxt("reducedelec_init.out", np.c_[energy["reducedelec"]])


for key in outcab:
    if "bl" not in key and 'err' not in key:
        np.savetxt(key+"d"+corr_name+"_"+mol_name+".dat", np.c_[time,outcab[key], outcab["err"+key]])
    if "bl" in key:
        for block in range(inputparam.nblocks):
            np.savetxt(key+str(block)+"d"+corr_name+"_"+mol_name+".dat", np.c_[time,outcab[key][block]])

print("combination complete")






