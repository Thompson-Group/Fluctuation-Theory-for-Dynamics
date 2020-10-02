#!/usr/bin/env python
"""
This is a general utility code that does a couple of things.
1) This grabs the activation energy decomposition from the fit_results files
2) This looks in the flucts.inp file for which energies exist
3) This calculates the elec component by subtracting ke and lj from tot
"""
import numpy as np
import sys

names = np.genfromtxt('flucts.inp', usecols=0, unpack=True,dtype=str)


print(names)

Eaint, Eatau, EaD = [],[],[]
errint, errtau, errD= [],[],[]
for name in names:
    if len(sys.argv) == 1 or sys.argv[1] == "c2":
        with open('fit_results_'+str(name)+'_'+str(name)+'_water_c2.dat') as f:
            lines = f.readlines()
            for line in lines:
                if "Ea Tau 2:" in line:
                    Eatau.append(float(line.split()[3]))
                    errtau.append(float(line.split()[5]))
                if "<EaTau2>" in line:
                    Eaint.append(float(line.split()[2]))
                    errint.append(float(line.split()[4]))
    else:
        Eatau.append(0.0)
        errtau.append(0.0)
        Eaint.append(0.0)
        errint.append(0.0)
    if len(sys.argv) == 1 or sys.argv[1] == "msd":
        with open('fit_results_'+str(name)+'_'+str(name)+'_water_msd.dat') as f:
            lines = f.readlines()
            for line in lines:
                if "Ea:" in line:
                    EaD.append(float(line.split()[1]))
                    errD.append(float(line.split()[3]))
    else:
        EaD.append(0.0)
        errD.append(0.0)

if "ewald" in names:
    if "coul" in names:
        eind = names.tolist().index("e")
        kind = names.tolist().index("ke")
        ljind = names.tolist().index("lj")
        names = names.tolist()
        names.append("elec")
        EaD.append(EaD[eind]-EaD[kind]-EaD[ljind])
        errD.append(np.sqrt(errD[eind]**2. + errD[kind]**2. + errD[ljind]**2.))
        Eaint.append(Eaint[eind]-Eaint[kind]-Eaint[ljind])
        errint.append(np.sqrt(errint[eind]**2. + errint[kind]**2. + errint[ljind]**2.))
        Eatau.append(Eatau[eind]-Eatau[kind]-Eatau[ljind])
        errtau.append(np.sqrt(errtau[eind]**2. + errtau[kind]**2. + errtau[ljind]**2.))


np.savetxt('decomp_data.dat', np.c_[names,EaD, errD, Eatau, errtau, Eaint, errint], fmt="%s")

   
