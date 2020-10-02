#!/usr/bin/env python
import numpy as np
import sys


def dist_components(rA,rB,q,rmax_mask,op):
    """
    This is a pretty complicated function to do a simple thing.
    This code takes two vectors of size(m,3) and size(n,3)
    Then it uses numpy broadcasting to calculate ALL the pairwise distances and outputs them in a matrix, sized (m,n,3)
    Then I use it to calculate the distances, and the vectors.
    (described here: https://stackoverflow.com/questions/60039982/numpy-python-vectorize-distance-function-to-calculate-pairwise-distance-of-2-ma/60040269#60040269)
    """
    vecdr = rA[:,np.newaxis,:]-rB[np.newaxis,:,:]
    vecdr = vecdr - np.multiply(L,np.round(np.divide(vecdr,L)))
    dr = np.linalg.norm(vecdr,axis=-1)
    if op == 0:
        return (np.abs(dr)<rmax)
    else:
        Efield = np.sum(coul(dr,vecdr,q,rmax_mask),axis=1)
        edr = np.divide(vecdr,dr[:,:,np.newaxis],where=dr[:,:,np.newaxis]!=0,out=np.zeros_like(vecdr))
        return edr, dr, Efield


def coul(dr,vdr, q,rmax_mask):
    conv = 1.0
    vdr=np.multiply(mask,vdr)
    num=np.multiply(vdr,q)
    tmpdr = dr*(rmax_mask)
    return conv*np.divide(num,np.power(tmpdr[:,:,np.newaxis],3),where=tmpdr[:,:,np.newaxis]!=0.0,out=np.zeros_like(num))

def E_Field(EOO,E1O,E2O,E12,E21,E11,E22,e1,e2):
    """
    Calculates Coulomb's law on the atoms
    """
    E1=np.add(np.add(E1O,E12),E11)
    E2=np.add(np.add(E2O,E21),E22)
    E1=np.multiply(E1,angperausq)
    E2=np.multiply(E2,angperausq)
    e1x,e1y,e1z=np.diagonal(np.swapaxes(e1,0,2)[0]),np.diagonal(np.swapaxes(e1,0,2)[1]),np.diagonal(np.swapaxes(e1,0,2)[2])
    e2x,e2y,e2z=np.diagonal(np.swapaxes(e2,0,2)[0]),np.diagonal(np.swapaxes(e2,0,2)[1]),np.diagonal(np.swapaxes(e2,0,2)[2])
    e1 = np.swapaxes(np.array((e1x,e1y,e1z)),0,1)
    e2 = np.swapaxes(np.array((e2x,e2y,e2z)),0,1)
    Efield1 = np.sum(np.multiply(e1,E1),axis=1)
    Efield2 = np.sum(np.multiply(e2,E2),axis=1)
    Efield = np.append(Efield1,Efield2)
    return Efield

def find_M(rO,r1,r2):
    """
    Subroutine to calculate the location of the point charge in tip4p2005
    Note that alpha is defined at the bottom of this code (and is 0.5*0.1546/(cos(0.5*ang)*blen))
    """
    tmpM = np.add(np.subtract(r1,rO),np.subtract(r2,rO))
    tmpM = np.multiply(alpha,tmpM)
    rM = np.add(rO,tmpM)
    return rM


def read_frames(f,t):
    """
    Reads the frames within filename.
    """
    # Defines output arrays
    rO,r1,r2 = [],[],[]
    eOO,eO1,eO1,dOO,dHO1,dHO2=[],[],[],[],[],[]
    count,index = 0,0
    while True:
        line = f.readline()
        if len(line.split()) == 1 and count > 0:
            if t != -1:
                rO,r1,r2 = np.array(rO),np.array(r1),np.array(r2)
                if tip4pflag==1: rO = find_M(rO,r1,r2) # If tip4p, replaces rO distance with rM distance
                # Calculate rmask
                o1mask=dist_components(r1,rO,qO,0,0)
                o2mask=dist_components(r2,rO,qO,0,0)
                rmax_mask=o1mask|o2mask
                tmpeOO,tmprOO,tmpEOO=dist_components(rO,rO,qO,rmax_mask,1)
                tmpe1O,tmpr1O,tmpE1O=dist_components(r1,rO,qO,rmax_mask,1) 
                tmpe2O,tmpr2O,tmpE2O=dist_components(r2,rO,qO,rmax_mask,1)
                tmpe12,tmpr12,tmpE12=dist_components(r1,r2,qH,rmax_mask,1) 
                tmpe21,tmpr21,tmpE21=dist_components(r2,r1,qH,rmax_mask,1) 
                tmpe11,tmpr11,tmpE11=dist_components(r1,r1,qH,rmax_mask,1) 
                tmpe22,tmpr22,tmpE22=dist_components(r2,r2,qH,rmax_mask,1)
                # These arrays take the shape of (mol1,mol2,3)
                dOO,dHO1,dHO2,dHH = tmprOO,tmpr1O,tmpr2O,tmpr12
                eO1,eO2   =  tmpe1O, tmpe2O
                Efield = E_Field(tmpEOO,tmpE1O,tmpE2O,tmpE12,tmpE21,tmpE11,tmpE22,eO1,eO2)
                return Efield
            else:
                return len(rO)
        if len(line.split()) == 4:
            # Check if oxygen
            if count % 3 == 0:
                # Read oxygen atoms and wrap boundary conditions
                rO.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            elif count % 3 == 1:
                # Read in hydrogens and wrap boundary conditions
                r1x = float(line.split()[1])
                r1x = r1x + L*round((rO[index][0]-r1x)/L)
                r1y = float(line.split()[2])
                r1y = r1y + L*round((rO[index][1]-r1y)/L)
                r1z = float(line.split()[3])
                r1z = r1z + L*round((rO[index][2]-r1z)/L)
                r1.append([r1x, r1y, r1z])
            elif count % 3 == 2:
                # Read in hydrogens and wrap boundary conditions
                r2x = float(line.split()[1])
                r2x = r2x + L*round((rO[index][0]-r2x)/L)
                r2y = float(line.split()[2])
                r2y = r2y + L*round((rO[index][1]-r2y)/L)
                r2z = float(line.split()[3])
                r2z = r2z + L*round((rO[index][2]-r2z)/L)
                r2.append([r2x, r2y, r2z])
                index += 1
            else:
                print("Incorrect count during read")
            count += 1
    print("read complete")
    return 



def calc_w(Efield):
    """
    Calculates the frequencies from the electric field
    """
    term2=np.multiply(c1,Efield)
    term3=np.multiply(c2,np.power(Efield,2))
    w=np.add(np.add(c0,term2),term3)
    return w 

def calc_x10(w):
    """
    This calculates the taylor expansion <1|x|0> of w
    """
    term2 = np.multiply(xc2,w)
    return np.subtract(xc1,term2)

def calc_mu_prime(Efield):
    term2 = np.multiply(mc2,Efield)
    term3 = 0.0
    if mc3 != 0.0: term3 = np.multiply(mc3,np.power(Efield,2))
    return np.add(np.add(mc1,term2),term3)

def calc_transition_dipole(Efield,w):
    """
    Calculates the transition dipole moment from the electric field and frequency
    """
    x10 = calc_x10(w)
    mu_prime_rat = calc_mu_prime(Efield)
    mu10 = np.multiply(x10,mu_prime_rat)
    return mu10



import time

tread = time.time()

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-avfreq', default=0.0, type=float, help='This is the average frequency in wavenumbers, if you are recalculating, otherwise defaults to 0')
parser.add_argument('-timecut', default=10000, type=int, help='Number of steps to cut calculation after')
args   = parser.parse_args()

avfreq = args.avfreq
timecut = args.timecut


# Read Corr Calc Input File
inpfile='corr_calc.in'
inp = open(inpfile, 'r')
lines = []
for line in inp:
    lines.append(line)
nfile  = str(lines[1]).strip()
ntimes = int(lines[3].split()[0])
dt     =float(lines[3].split()[1])
volume = float(lines[5])
mol_name = str(lines[7]).strip()
qO = float(lines[9])
inp.close()



tip4pflag = 0 # 0 - not tip4p, 1 tip4p
alpha=0.0 # this is for finding the location of the point charge in TIP4P/2005
rmax=0.0
c0,c1,c2=0.0, 0.0, 0.0
if qO == -0.84760:
    print("Choosing SPC/E for water model")
    qO = -0.84760
    rmax=7.831
    # w coeffs
    c0=3761.6
    c1=-5060.4
    c2=-86225.0
    # x_10 coeffs 
    xc1=0.71116 
    xc2=0.927E-5
    # mu_10 coeffs
    mc1=0.71116
    mc2=75.591
    mc3=0.0
elif qO == -1.1128:
    print("Recognized as TIP4P/2005 - Note This does not have a map currently")
    tip4pflag = 1
    hoh_ang = 104.52
    blength = 0.9572
    alpha = 0.1546/(2*np.cos(0.5*hoh_ang*np.pi/180.0)*blength)
    rmax=7.831
    # w coeffs
    c0=3760.2
    c1=-3541.7
    c2=-152677
    # x_10 coeffs
    xc1=0.19285
    xc2=1.7261E-5
    # mu_10 coeffs
    mc1=0.1646
    mc2=11.39
    mc3=63.41
elif qO == -1.040:
    print("Recognized as TIP4P")
    tip4pflag = 1
    hoh_ang = 104.52
    blength = 0.9572
    alpha = 0.15/(2*np.cos(0.5*hoh_ang*np.pi/180.0)*blength)
    rmax=7.831
    # w coeffs
    c0=3760.2
    c1=-3541.7
    c2=-152677
    # x_10 coeffs
    xc1=0.19285
    xc2=1.7261E-5
    # mu_10 coeffs
    mc1=0.1646
    mc2=11.39
    mc3=63.41

print(c0,c1,c2)

qH=-qO/2
L = volume**(1/3.)
angperau=0.52917721092
angperausq  = np.power(angperau,2)

# Trajectory File Name
filename = 'traj_'+str(nfile)+'_'+str(mol_name)+'.xyz'

# Open Traj, read number of molecules, then close and reopen
f=open(filename,'r')
nmols = int(read_frames(f,-1))
f.close()
f=open(filename,'r')

# Header Info
# Sets the mask for making the like indices zero
vr = np.array((np.linspace(1,nmols,num=nmols),np.linspace(1,nmols,num=nmols),np.linspace(1,nmols,num=nmols)))
vr=np.swapaxes(vr,0,1)
mask=vr[:,np.newaxis,:]-vr[np.newaxis,:,:]
mask=1*(mask!=0)
# Sets up the zero field

if ntimes > timecut:
    ntimes = timecut

Efield=np.zeros((ntimes,2*nmols))
# calculates for t=0
# Loop over times
for t in range(ntimes):
    if t%100 == 0: print("STEP:",t,flush=True)
    Efield[t]= read_frames(f,t)
# Converts to frequencies


w = calc_w(Efield)
w_avg = avfreq
dw = np.subtract(w,w_avg)
w2avg = np.average(w)
dw2 = np.subtract(w,w2avg)

print("Finalizing Calculations",flush=True)

# Calculates transition dipole
initial_w = w[0]
initial_mu10 = calc_transition_dipole(Efield[0],w[0])

# Calculates the time correlation function
corr=np.average(np.multiply(dw[0],dw),axis=1)
corr2=np.average(np.multiply(dw2[0],dw2),axis=1)


# Outputs the final data
import pickle

pickle.dump(corr,open('spec_diff_'+nfile+'_'+mol_name+'.pckl','wb'))
pickle.dump(corr2,open('spec_diff2_'+nfile+'_'+mol_name+'.pckl','wb'))
pickle.dump(initial_w,open('spec_dist_'+nfile+'_'+mol_name+'.pckl','wb'))
pickle.dump(initial_mu10, open('spec_mu10_'+nfile+'_'+mol_name+'.pckl','wb'))

print("End time %2.5f" % (time.time()-tread))
