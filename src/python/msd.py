#!/usr/bin/env python
"""
This is a modified version of my msd code to work with the fluctuation method. Note - it calculates multiple origins and thus can't be used for a decomposition; however, it is powerful for calculating the derivative of the diffusion coefficient with more accuracy.

Copyright Zeke Piskulich 2020
"""
import numpy as np
import pickle
import argparse,gzip

def choose_mol_read(i):
    def mol_read_unwrapped(N,AperMol, Mcomp,f):
        mtot = np.sum(Mcomp)
        rcom = []
        rcomsys = np.zeros(3)
        for i in range(N):
            rtmp = np.zeros((AperMol,3))
            for j in range(AperMol):
                line = f.readline().strip().split()
                r    = np.array(line[1:],dtype=float)
                rtmp[j] = Mcomp[j]/mtot*r
            rcom.append(np.sum(rtmp,axis=0))
            
        return rcom
    def mol_read_wrapped(N,AperMol, Mcomp,f):
        mtot = np.sum(Mcomp)
        rcom = []
        for i in range(N):
            rtmp = np.zeros((AperMol,3))
            r0 = np.zeros(3)
            for j in range(AperMol):
                line = f.readline().strip().split()
                r    = np.array(line[1:],dtype=float)
                if j != 0:
                    r = r+L*np.round((r0-r)/L)
                else:
                    r0 = r
                rtmp[j] = Mcomp[j]/mtot*r
            rcom.append(np.sum(rtmp,axis=0))
        return rcom
    if i == 0:
        return mol_read_unwrapped
    else:
        return mol_read_wrapped

def read_molfile(mol):
    mfile = open("../../"+mol+".txt", 'r')
    atoms = int(mfile.readline().strip().split()[0])
    mfile.readline()
    mfile.readline()
    m = []
    for i in range(atoms):
        m.append(float(mfile.readline().strip().split()[0]))
    return atoms, m
    

def read_traj(filename):
    if ".gz" in filename:
        f = gzip.open(filename,'r')
    else:
        f = open(filename,'r')
    rcom = {}
    mol_read = choose_mol_read(unwrap)
    nlines = 0
    for component in components:
        nlines += N[component]*AperMol[component]
        rcom[component] = np.zeros((nframes,N[component],3))
    nlines += skip
    for frame in range(start_conf):
        if frame%1000 == 0: print("skipping frame %d" % frame)
        for line in range(nlines):
            f.readline()
    rcom_sys = np.zeros((nframes,3))
    Mtotsys = 0.0
    for component in components:
        Mtotsys += np.sum(M[component])*N[component]
    for frame in range(nframes):
        if frame%1000 == 0: print(frame)
        for i in range(skip):
            f.readline()
        for component in components:
            rcom[component][frame]=mol_read(N[component],AperMol[component],M[component],f)
            rcom_sys[frame] += np.sum(rcom[component][frame]*np.sum(M[component]),axis=0)/Mtotsys
        rcom[component][frame] -= rcom_sys[frame]
    return rcom

def choose_calc_msd(i):
    def calc_msd_unwrapped(r,n,mol):
        last_t0 = nframes-corr_len
        msd = []
        nt0 = int(last_t0/or_split)
        
        for t0 in range(0,last_t0,or_split):
            tmp = np.zeros(n)
            tmpmsd = []
            for j in range(t0,t0+corr_len):
                drsq = np.sum(np.power(r[j]-r[t0],2),axis=-1)
                tmp = drsq
                tmpmsd.append(np.average(tmp))
            msd.append(tmpmsd)
        x = np.arange(nframes)
        time = np.multiply(x,timestep)[:corr_len]
        if write_pckl==1: pickle.dump(msd,open('msd.pckl','wb'),protocol=4)
        np.savetxt("msd_"+mol+"_out.dat",np.c_[time,np.average(msd,axis=0)])
        return time, np.average(msd,axis=0)
    def calc_msd_wrapped(r,n,mol):
        last_t0 = nframes-corr_len
        msd = []
        for t0 in range(0,last_t0,or_split):
            shift,tmp = np.zeros((n,3)),np.zeros(n)
            tmpmsd = []
            for j in range(t0,t0+corr_len):
                if j!=t0: shift -= L*np.round((r[j]-r[j-1])/L)
                drsq = np.sum(np.power(r[j]-r[t0] + shift,2),axis=-1)
                tmp = drsq
                tmpmsd.append(np.average(tmp))
            msd.append(tmpmsd)
        x = np.arange(nframes)
        time = np.multiply(x,timestep)[:corr_len]
        if write_pckl==1: pickle.dump(msd,open('msd.pckl','wb'),protocol=4)
        np.savetxt("msd_"+mol+"_out.dat",np.c_[time,np.average(msd,axis=0)])
        return time, np.average(msd,axis=0)
    if i==0: 
        return calc_msd_unwrapped
    else:
        return calc_msd_wrapped

def fit_msd(time,msd):
    def linear(x,m,b):
        return m*x + b
    from scipy.optimize import curve_fit
    def calc_D(m):
        return m/0.6
    cut = int(len(msd)*0.75)
    popt,pcov = curve_fit(linear, time[cut:], msd[cut:])
    m,b = popt
    return calc_D(m) 

         

parser = argparse.ArgumentParser()
parser.add_argument('-fname', default="traj.xyz", type=str, help="The file name for the trajectory")
parser.add_argument('-comp', default=[], action='append', help='This is the name of each component (repeatable)')
parser.add_argument('-N',  default=[], action='append', help="The number of each component (repeatable)")
parser.add_argument('-corr_len', default=1000, type=int, help="Correlation length")
parser.add_argument('-or_split', default=20, type=int, help="Separation of origins")
parser.add_argument('-timestep', default = 0.1, type=float, help="Timestep (in ps)")
parser.add_argument('-skip', default=2, type=int, help="Lines to skip at beginning of each frame")
parser.add_argument('-nframes', default=50000, type=int, help="Number of total frames")
parser.add_argument('-Lfile', default="NL.avg", type=str, help="File with box length")
parser.add_argument('-op', default=0, type=int, help="Options")
parser.add_argument('-unwrap', default=1, type=int, help="[1] if need to unwrap, [0] if already unwrapped")
parser.add_argument('-start_conf',default=0, type=int, help="Starting configuration")
parser.add_argument('-write_pckl', default=0, type=int, help="Write a pckl file with configs, [1] yes, [0] no")
args = parser.parse_args()

filename   = args.fname
components = args.comp
n          = np.array(args.N,dtype=int)
corr_len   = args.corr_len
or_split   = args.or_split
timestep   = args.timestep
skip       = args.skip
nframes    = args.nframes
Lfname     = args.Lfile
op         = args.op
unwrap     = args.unwrap
start_conf = args.start_conf
write_pckl = args.write_pckl

global f

filename = "traj_"

L=0.0
num = 0
with open("corr_calc.in",'r') as lf:
    lf.readline()
    num=int(lf.readline().strip().split()[0])
    lf.readline()
    lf.readline()
    lf.readline()
    L = float(lf.readline().strip().split()[0])**(1/3.)
print(L)

filename = "traj_"+str(num)+"_"+components[0]+".xyz"
print(filename)

AperMol, M, N  = {},{},{}
count = 0
for component in components:
    AperMol[component], M[component] = read_molfile(component)
    N[component] = n[count]
    count +=1

RCOM={}
if op == 0:
    print("Reading Trajectory")
    RCOM = read_traj(filename)
    if write_pckl ==1: pickle.dump(RCOM,open('rcom.pckl','wb'),protocol=4)
    print("Center of Mass Coordinates Recorded")
else:
    print("Reading COM")
    RCOM = pickle.load(open('rcom.pckl','rb'))
    print("COM is Read")

calc_msd = choose_calc_msd(unwrap)

msd, D = {},{}
time = []
for component in components:
    time, msd[component]=calc_msd(RCOM[component],N[component],component)
    D[component]=fit_msd(time,msd[component])
    print("Diffusion coefficient  of %s is %s x 10^5 cm^2/s" % (component,D[component]))
    pickle.dump(D[component],open('D_'+str(num)+'_'+component+'.pckl','wb'))
    pickle.dump(msd[component],open('msd_py_'+str(num)+'_'+component+'.pckl','wb'))
           
