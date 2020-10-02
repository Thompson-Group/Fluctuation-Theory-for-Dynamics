#!/usr/bin/env python
#This is a code aimed at reselecting the velcoties from the appropriate temperature distribution.
"""
This is used in the Ion Pairing code calculation. This takes the current velocities and then creates a 
velocity in the axial direction i.e. to undo rigid bonds.
"""
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import sys
from shutil import copyfile


def pick_mb(mean, std, npick):
    # Python Function to Pull Samples from a Distribution
    samples=npr.normal(mean, std, npick)
    return samples[0]

def plot_hist(data, bins,mean,std):
    count, bins, ignored = plt.hist(data, bins, normed=True)
    plt.plot(bins, 1/(std * np.sqrt(2 * np.pi)) * np.exp( - (bins-mean)**2 / (2 * std**2) ))
    plt.show()
    return None

def read_restart(filename,skeyword, ekeyword):
    with open(filename, 'r') as f_restart:
        count = 0
        first_atom=0
        lines=[]
        for line in f_restart:
            count += 1
            if skeyword in line and count > 50:
                for line in f_restart:
                    if ekeyword in line:
                        break
                    lines.append(line.strip())
    return lines

def search_string(lines, key):
    count = 0
    for line in lines:
        count += 1
        if key in line:
            break
    return count

def searchfile(filename):
    searchfile = open(filename, 'r')
    count = 0
    for line in searchfile:
        count += 1
        if "&END VELOCITY" in line and count > 50:
            break
    searchfile.close()
    return count
    
            

def grab_vals(lines, keyval):
    vals = lines[keyval].split()
    vals = [float(i) for i in vals]
    return vals

def calc_com_vel(v,i):
    vcom = (float(v[0][i])*mLi+float(v[1][i])*mF)/M
    return vcom

def calc_rel_vel(v,i):
    vrel = float(v[0][i])-float(v[1][i])
    return vrel

def calc_delv(vrel,vrelo):
    delvF=(vrel-vrelo)*mLi/M
    delvLi=-(vrel-vrelo)*mF/M
    return delvLi, delvF

def write_restart(filename, vo, vn, line1):
    count = 0
    line2 = line1+1
    with open(filename, 'r') as f:
        out = open('tmp', 'w')
        for line in f:
            count += 1
            if count == line1:
                print("test")
                out.write("           %2.16e   %.16e    %.16e\n" % (v[0][0], v[0][1], v[0][2]))
            elif count == line2:
                out.write("           %2.16e   %.16e    %.16e\n" % (v[1][0], v[1][1], v[1][2]))
            else:
                out.write(line)
        out.close()
    return None
        


#samples=pick_mb(0.0, 1.0, 2000)
#plot_hist(samples, 200, 0.0, 1.0)

# Array Definitions
v=[]
vo=[[0,0,0],[0,0,0]]
vcom=[]

# Constants
kb=1.38064852E-23 # J/K
bohr_in_m = 5.2918E-11  # m
au_time_in_fs = 0.0242  # fs
fs_in_s = 1E-15 # s
conv_vel_to_cp2k = (1.0/bohr_in_m)*fs_in_s*au_time_in_fs # bohr/au_time
Na=6.0221409E23 # molecs
g_in_kg = 1E-3 # kg


# Command Line Input
if(len(sys.argv) != 5):
    print("Usage vel_reselect.py mLi mF T")
    exit(0)
mLi=float(sys.argv[1])
mF=float(sys.argv[2])
T=float(sys.argv[3])
filename=str(sys.argv[4])

# Calculated Quantities
M = mLi + mF
mu = mLi*mF/M

# Parse File into array
coords=read_restart(filename, "&COORD", "&END")
vels=read_restart(filename, "&VELOCITY", "&END")
val = search_string(coords, "Li")
v.append(grab_vals(vels,val-1))
v.append(grab_vals(vels,val))

# Calculate Current Com Vel
for i in range(0,3):
    vcom.append(calc_com_vel(v, i))
# Pull a velocity from the Distribution
mean=0.0
sigma=np.sqrt(kb*T/mu*(Na)*(1/g_in_kg))
for i in range(0,3):
    vrelo = calc_rel_vel(v,i)
    vrel = pick_mb(mean, sigma, 1)
    vrel *= conv_vel_to_cp2k
    delv= calc_delv(vrel,vrelo)
    delvLi=delv[0]
    delvF=delv[1]
    vold = v[0][i]
    vo[0][i]=vold
    vold = v[1][i]
    vo[1][i]=vold
    v[0][i] += delvLi
    v[1][i] += delvF
    print("delvF %s delvLi %s" % (delvF, delvLi))
    print("Com %s Orig %s" % (calc_com_vel(v,i),vcom[i]))
    print("Rel %s" % calc_rel_vel(v,i))
line=searchfile(filename) - 2
print(line)
write_restart(filename, vo, v,line)
copyfile("tmp",filename)


