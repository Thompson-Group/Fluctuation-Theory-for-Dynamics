#!/usr/bin/env python
""" 
I realize this isn't a fortran program - but it was much simpler to write.
"""
import numpy as np
import time,sys
import argparse
from scipy.spatial.distance import pdist, squareform,cdist


# Input Parameters

parser=argparse.ArgumentParser()
parser.add_argument('-late_hbond', default=2000,type=int,help="At what timestep should jumps be split into early or late.")
args = parser.parse_args()

late_hbond = args.late_hbond

# SPC/E Water Model

qh = 0.42380 # e
qo = -0.84760 # e
eps = 0.15535 # Kcal/mol
sig = 3.16600 # Angstroms
doh = 1.0 # Angstroms

def parse_log(filename):
    """
    This subroutine reads the lammps log file and generates the trajectories.
    """
    count,start,end,empty=0,0,0,0
    lines=[]
    header=""
    with open(filename,'r') as f:
        lines=f.readlines()
        for line in lines:
            if "Step Time" in line:
                start=count
                header=line.strip().split()
            if "Loop time" in line:
                end=count
            if end != 0 and line == '\n':
                empty+=1
            count += 1

    startskip = start + 1
    endskip   = len(lines)-end-empty
    data = np.genfromtxt(filename, skip_header=startskip, skip_footer=endskip,comments=None,unpack=True)
    edata={}
    for i in range(len(header)):
        edata[header[i]]=data[i]
    return edata

def LJ(dr,eps,sig):
    """
    This function calculates the lennard jones potential energy for a particular epsilon and sigma
    """
    r6 = (sig/dr)**6.
    Ulj = 4*eps*(r6**2.-r6)
    return Ulj

def Coul(dr,q1,q2):
    """
    This function calculates the real-space electrostatic energy
    """
    C=332.06371
    Ucoul = C*q1*q2/dr
    return Ucoul

def calc_lj(mol1,mol2,eps,sig):
    """
    This calls the calc lj script
    """
    #dr = min_dist(rO[mol1],rO[mol2])[0]
    dr =1.0
    return LJ(dr,eps,sig)

def calc_coul(mol1,mol2,qo,qh):
    """
    This calls the electrostatic calculation script
    """
    terms = []
    # Note - there are more terms because there are many sites which have interactions! 
    # 2 water molecules, 3 atoms each - 9 total terms
    # The donor and acceptor are tagged     
    """D  A
    terms.append([min_dist(rO[mol1],rO[mol2])[0],qo,qo]) # Ox Ox
    terms.append([min_dist(rO[mol1],r1[mol2])[0],qo,qh]) # Ox H1
    terms.append([min_dist(rO[mol1],r2[mol2])[0],qo,qh]) # Ox H2
    terms.append([min_dist(r1[mol1],rO[mol2])[0],qh,qo]) # H1 Ox
    terms.append([min_dist(r2[mol1],rO[mol2])[0],qh,qo]) # H2 Ox
    terms.append([min_dist(r1[mol1],r1[mol2])[0],qh,qh]) # H1 H1
    terms.append([min_dist(r2[mol1],r1[mol2])[0],qh,qh]) # H2 H1
    terms.append([min_dist(r1[mol1],r2[mol2])[0],qh,qh]) # H1 H2
    terms.append([min_dist(r2[mol1],r2[mol2])[0],qh,qh]) # H2 H2"""
    Ucoul=0
    """
    # Sums up the terms
    for term in terms:
        Ucoul += Coul(term[0],term[1],term[2])
    """
    return Ucoul


def min_dist(a, b):
    """
    This calculates the minimum image distance between two 3D coordinates, stored 
    in a and b. This outputs the overall final distance as a float.
    """
    dr      = []
    shift   = []
    dsq     = 0.0
    for i in range(3):
        shift.append(L*round((a[i]-b[i])/L))
        dr.append(a[i]-b[i]-shift[i])
        dsq += dr[i]**2.
    drfinal = np.sqrt(dsq)
    return drfinal, dr


def calc_ang(eOO,eOH,dOO,dOH):
    """
    This calculates the angle between the Od->Oa vector and the Od->Hd vector
    This gives the angle of the hydrogen bond
    Note - output is in radians!
    """
    edote = np.dot(eOO,eOH)
    dHOO = np.arccos(edote)
    return dHOO


def dist_components(rA,rB):
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
    edr = np.divide(vecdr,dr[:,:,np.newaxis],where=dr[:,:,np.newaxis]!=0,out=np.zeros_like(vecdr))
    return edr, dr


def read_frames(f,t):
    """
    Reads the frames within filename.
    """
    # Defines output arrays
    rO,r1,r2 = [],[],[]
    eOO,eO1,eO1,dOO,dHO1,dHO2=[],[],[],[],[],[]
    count,nmols,index = 0,0,0
    while True:
        line = f.readline()
        if len(line.split()) == 1 and count > 0:
            rO,r1,r2 = np.array(rO),np.array(r1),np.array(r2)
            tmpeOO,tmprOO=dist_components(rO,rO)
            tmpe1O,tmpr1O=dist_components(r1,rO)
            tmpe2O,tmpr2O=dist_components(r2,rO)
            # These arrays take the shape of (time,mol1,mol2)
            dOO,dHO1,dHO2 = tmprOO,tmpr1O,tmpr2O
            eOO,eO1,eO2   = -tmpeOO, tmpe1O, tmpe2O # The minus sign is intentional - by default it gives the wrong direction - need to reverse it!
            return nmols, dOO, dHO1, dHO2, eOO, eO1, eO2
        if len(line.split()) == 4:
            # Check if oxygen
            if count % 3 == 0:
                # Read oxygen atoms and wrap boundary conditions
                rO.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
            elif count % 3 == 1:
                # Read in hydrogens and wrap boundary conditions
                r1x = float(line.split()[1])
                r1x = r1x - L*round((rO[index][0]-r1x)/L)
                r1y = float(line.split()[2])
                r1y = r1y - L*round((rO[index][1]-r1y)/L)
                r1z = float(line.split()[3])
                r1z = r1z - L*round((rO[index][2]-r1z)/L)
                r1.append([r1x, r1y, r1z])
            elif count % 3 == 2:
                # Read in hydrogens and wrap boundary conditions
                r2x = float(line.split()[1])
                r2x = r2x - L*round((rO[index][0]-r2x)/L)
                r2y = float(line.split()[2])
                r2y = r2y - L*round((rO[index][1]-r2y)/L)
                r2z = float(line.split()[3])
                r2z = r2z - L*round((rO[index][2]-r2z)/L)
                r2.append([r2x, r2y, r2z])
                index += 1
                nmols += 1
            else:
                print("Incorrect count during read")
            count += 1
    print("read complete")
    return nmols, dOO, dHO1, dHO2, eOO,eO1,eO2

def is_it_hbonded(mol1, mol2, dOO, dHO1, dHO2,eOO,eO1,eO2, t, orig):
    """
    Outputs whether it is hbonded or not.

    This works based on three criteria:
    1) is the OO distance < rOO_max
    2) is the HO distance < rOH_max
    3) is the HOO angle   < ang_max

    If yes to all three, outputs 1 or 2 depending on which hydrogen is involved in the hbond,
    otherwise, outputs 0.
    """
    hbonded = 0 # default not hbonded
    # Calculates OO Distance
    if dOO[t][mol1][mol2] == 0.0: print(dOO[mol1][mol2],mol1,mol2)
    if dOO[t][mol1][mol2] < rOO_max:
        # Calculates the distance of each OH on molecule 1 from the Oxygen
        if orig == 1 or orig == -1:
            if dHO1[t][mol1][mol2] < rHO_max:
                # Calculates the Angle of HOO from the donor and acceptor (first OH)
                dHOO1 = np.degrees(calc_ang(eOO[t][mol1][mol2],eO1[t][mol1][mol1],dOO[t][mol1][mol2], dHO1[t][mol1][mol1]))
                if dHOO1 < ang_max:
                    hbonded = 1
        if orig == 2 or orig == -1:
            if dHO2[t][mol1][mol2] < rHO_max:
                # Calculates the Angle of HOO from the donor and acceptor (second OH)
                dHOO2 = np.degrees(calc_ang(eOO[t][mol1][mol2],eO2[t][mol1][mol1],dOO[t][mol1][mol2], dHO2[t][mol1][mol1]))
                if dHOO2 < ang_max:
                    hbonded = 2
    return hbonded

def calc_cn(e_o, e_t):
    """
    This calculates the value of the reorientational correlation functions (orders = 1-3)
    These are described in successive orders of the legendre polynomials
    P1 = x
    P2 = 0.5(3x^2 - 1)
    P3 = 0.5(5x^3 - 3x)
    
    Output: A tuple of the results of all 3.
    """
    edote = np.dot(e_o, e_t)
    c1 = edote
    c2 = 0.5*(3*edote**2. - 1)
    c3 = 0.5*(5*edote**3-3*edote)
    return c1,c2,c3

def calc_jumpang(OH,eOO,t):
    """
    Calculates the jump angle after a switch
    
    This is defined as the angle between the vectors
    Od->Oa0 (original acceptor angle at time of jump)
    Od->Oan (new acceptor angle at the time of the jump)
    Output: the angle of the jump in radians.
    """
    mol1, molo, mol2 = OH[0],OH[1],OH[6]
    e_old=eOO[t][mol1][molo]
    e_new=eOO[t][mol1][mol2]
    theta = np.arccos(np.dot(e_old, e_new))
    return theta


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
inp.close()

L = volume**(1/3.)

# Definition of the HBond Criteria
rOO_max = 3.1
rHO_max = 2.0
ang_max = 20.

# Trajectory File Name
filename = 'traj_'+str(nfile)+'_'+str(mol_name)+'.xyz'

# Read the frames
print('Reading frames')
tread=time.time()
f=open(filename,'r')
nmols, dOOo, dHO1o, dHO2o, eOOo,eO1o,eO2o= read_frames(f,0)
dOO=np.array([dOOo,dOOo])
dHO1=np.array([dHO1o,dHO1o])
dHO2=np.array([dHO2o,dHO2o])
eOO=np.array([eOOo,eOOo])
eO1=np.array([eO1o,eO1o])
eO2=np.array([eO2o,eO2o])
print('Initial Frame Reading Complete: %2.5f seconds' % (time.time()-tread))
# Determine initial HBONDS (time 0)
"""
Form of the OHs vector:
    This has multiple elements that are important
    The first index provides the row of interest, corresponding to a particular OH
    The second index provides the column of interest which is separated into a few different parts
    moloh mola,old on? moldonor1 moldonor2 mola,new
"""
OHs = []

print('Calculating initial hbonds')
for mol1 in range(nmols):
    for mol2 in range(nmols):
        if mol1 != mol2:
            hbnd = is_it_hbonded(mol1, mol2, dOO, dHO1,dHO2,eOO,eO1,eO2,0,-1) # hbnd has options 0 (not hbonded) 1 (first oh hbonded) and 2 (second oh hbonded)
            if hbnd != 0:
                OHs.append([mol1,mol2,hbnd])


# This section checks for donors (up to a maximum of 3)
for OH1 in range(len(OHs)):
    count = 0
    if count < 2: # Maximum of 3 donors planned
        for OH2 in range(len(OHs)):
            if OHs[OH2][0] == OHs[OH1][1]:
                OHs[OH2].append(OHs[OH1][0])
                count+=1
# This section fills out the donors with -1's where there are no donors
for OH1 in range(len(OHs)):
    if len(OHs[OH1]) == 3:
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
    elif len(OHs[OH1]) == 4:
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
    elif len(OHs[OH1]) == 5:
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
    elif len(OHs[OH1]) == 6:
        OHs[OH1].append(-1)


print("There are %s OHs" % len(OHs))
print("%d hydrogen bonds" % np.sum(np.array(OHs).T[2]>0))

# Sets frequency to print progress.
nval = int(ntimes/10.)

crp = np.zeros((len(OHs),ntimes))                                                                       # This is the jump correlation function (side-side TCF)
c1,c2,c3  = np.zeros((len(OHs), ntimes)),np.zeros((len(OHs), ntimes)),np.zeros((len(OHs), ntimes))      # These are the frame reorientation correlation functions Od-Oh
ch1,ch2,ch3 = np.zeros((len(OHs), ntimes)),np.zeros((len(OHs), ntimes)),np.zeros((len(OHs), ntimes))    # These are the frame reorientation correlation functions OHd
norm_t = np.zeros(ntimes)                                                                               # This is the hydrogen bond extinction correlation function
theta,jindex = [],[]                                                                                    # This is the jump angle distribution - note NOT A TCF

c1_1, c1_2, c1_3, c1_4, c1_5 = np.zeros(ntimes), np.zeros(ntimes), np.zeros(ntimes), np.zeros(ntimes), np.zeros(ntimes)

# Setup Data Structure for Separating Frame Contributions
is_hbond_late = np.zeros(len(OHs))  # fills with 1s or 0s, default is not




# Loop over time, checks if hbonded at each time
"""
This section calls the loop over times.

Note, the loop goes over times, then loops over OHs, and then loops over every other molecule.
There are separate correlation functions calculated for the frame reorientation of the OO and the OH.
"""

t0=time.time()
end=0

# Defines the oo distance calculator
oo_dist = []

for n in range(1,ntimes): # Loop Over Times
    nmols, dOO[1], dHO1[1], dHO2[1], eOO[1],eO1[1],eO2[1]= read_frames(f,0)   
    if n%nval == 0: # Outputs information to the screen, what step is reached, how much time has elapsed, and the # of hbonds
        t1 = time.time()
        print("Step Reached: %2.2f (ps), Av Step Time %s seconds" % (n*dt,(t1-t0)/float(nval)))
        print("%d hydrogen bonds" % np.sum(np.array(OHs).T[2]>0))
        t0 = time.time()
    for OH in range(len(OHs)): # Loop over OHs.
        crp[OH][n] = crp[OH][n-1] # Sets the value of crp to the value at the previous step
        mol1 = OHs[OH][0]
        for mol2 in range(nmols):
            hbnd = 0
            if mol1 != mol2 and OHs[OH][2] != 0:
                hbnd = is_it_hbonded(mol1, mol2, dOO, dHO1, dHO2,eOO,eO1,eO2, 1,OHs[OH][2])
                if hbnd != 0:
                    if OHs[OH][1] == mol2: # Same Acceptor
                        crp[OH][n] += 0
                        # Calculates the oo c2
                        c1[OH][n],c2[OH][n],c3[OH][n]=calc_cn(eOO[0][mol1][mol2],eOO[1][mol1][mol2])
                        if OHs[OH][2] == 1:
                            ch1[OH][n],ch2[OH][n],ch3[OH][n]=calc_cn(eO1[0][mol1][mol1],eO1[1][mol1][mol1])
                        else:
                            ch1[OH][n],ch2[OH][n],ch3[OH][n]=calc_cn(eO2[0][mol1][mol1],eO2[1][mol1][mol1])
                        # Time dependent normalization
                        norm_t[n] += 1.0
                    else: # New Acceptor
                        if n >= late_hbond: is_hbond_late[OH] = 1.0   # This sets a boolean element to 1 if the hbond lasts a very long time.
                        crp[OH][n] += 1
                        c1[OH][n],c2[OH][n],c3[OH][n]=0.0,0.0,0.0
                        ch1[OH][n],ch2[OH][n],ch3[OH][n]=0.0,0.0,0.0
                        OHs[OH][2]=0
                        OHs[OH][6]=mol2
                        oo_dist.append(dOO[1][mol1][OHs[OH][1]])
                        theta.append(calc_jumpang(OHs[OH],eOO,1))
                        jindex.append(n)
                if hbnd == 0 and mol2 == OHs[OH][1]:
                    c1[OH][n],c2[OH][n],c3[OH][n]=calc_cn(eOO[0][mol1][mol2],eOO[1][mol1][mol2])
                    if OHs[OH][2] == 1:
                        ch1[OH][n],ch2[OH][n],ch3[OH][n]=calc_cn(eO1[0][mol1][mol1],eO1[1][mol1][mol1])
                    else:
                        ch1[OH][n],ch2[OH][n],ch3[OH][n]=calc_cn(eO2[0][mol1][mol1],eO2[1][mol1][mol1])
                    norm_t[n] += 1.0 
    if np.sum(np.array(OHs).T[2]>0) == 0:
        end = n+1
        break
    else: end = -1
f.close()
# This ends the time loop

# This is added to make sure that the correlations are calculated even past 
# the point at which there are no longer hydrogen bonds (just sets the remaining elements to 1)
if end != -1:
    for n in range(end,ntimes):
        for OH in range(len(OHs)):
            crp[OH][n] += 1
# This handles the case where there is no new acceptor 
for OH in range(len(OHs)):
    if len(OHs[OH]) != 7:
        OHs[OH].append(-1) 

# This calculates the jump correlations based on whether they are hydrogen bonded or not.
c1_early,c1_late=c1[is_hbond_late==0.0],c1[is_hbond_late==1.0]
num_early,num_late = np.shape(c1_early)[0],np.shape(c1_late)[0]
early_norm = norm_t - num_late
EC1 = np.divide(np.sum(c1_early,axis=0),early_norm, where=early_norm!=0,out=np.zeros_like(np.sum(c1_early,axis=0)))
LC1 = []
late_norm = norm_t - num_early
if num_late != 0:
    LC1 = np.divide(np.sum(c1_late,axis=0),late_norm, where=late_norm!=0, out=np.zeros_like(np.sum(c1_late,axis=0)))
else:
    LC1 = np.zeros(ntimes)

EC1[0],LC1[0] = 0.0,0.0

print("Jiggying up final calculations")

# This section averages all the frame contributions by dividing by the time dependent normalization.
C1 = np.divide(np.sum(c1,axis=0),norm_t,out=np.zeros_like(np.sum(c1,axis=0)),where=norm_t!=0)
C2 = np.divide(np.sum(c2,axis=0),norm_t,out=np.zeros_like(np.sum(c2,axis=0)),where=norm_t!=0)
C3 = np.divide(np.sum(c3,axis=0),norm_t,out=np.zeros_like(np.sum(c3,axis=0)),where=norm_t!=0)
CH1 = np.divide(np.sum(ch1,axis=0),norm_t,out=np.zeros_like(np.sum(ch1,axis=0)),where=norm_t!=0)
CH2 = np.divide(np.sum(ch2,axis=0),norm_t,out=np.zeros_like(np.sum(ch2,axis=0)),where=norm_t!=0)
CH3 = np.divide(np.sum(ch3,axis=0),norm_t,out=np.zeros_like(np.sum(ch3,axis=0)),where=norm_t!=0)
# Sets the first element (which isn't otherwise calculated) to the correct initial value
C1[0],C2[0],C3[0]=1.0,1.0,1.0
CH1[0],CH2[0],CH3[0]=1.0,1.0,1.0
CRP = np.average(crp,axis=0)


# This section is for the special decompositions, currently is disabled.
UCRP={"LJAold":np.zeros((len(OHs),ntimes)), "LJAnew":np.zeros((len(OHs),ntimes)), "LJD":np.zeros((len(OHs),ntimes)),"CAold":np.zeros((len(OHs),ntimes)),"CAnew":np.zeros((len(OHs),ntimes)), "CD":np.zeros((len(OHs),ntimes))}
UC2={"LJAold":np.zeros((len(OHs),ntimes)), "LJAnew":np.zeros((len(OHs),ntimes)), "LJD":np.zeros((len(OHs),ntimes)),"CAold":np.zeros((len(OHs),ntimes)),"CAnew":np.zeros((len(OHs),ntimes)), "CD":np.zeros((len(OHs),ntimes))}

Ulj={"Aold":[],"Anew":[], "D":[]}
UCoul={"Aold":[],"Anew":[], "D":[]}

"""
for OH in range(len(OHs)):
    target, oldaccept, donor1, donor2, donor3, newaccept = OHs[OH][0],OHs[OH][1],OHs[OH][3],OHs[OH][4],OHs[OH][5], OHs[OH][6]
    # Old Acceptor
    Ulj["Aold"].append(calc_lj(target,oldaccept,eps,sig))
    UCRP["LJAold"][OH]=np.multiply(Ulj["Aold"][OH],crp[OH])
    UC2["LJAold"][OH]=np.multiply(Ulj["Aold"][OH],c2[OH])
    UCoul["Aold"].append(calc_coul(target,oldaccept,qo,qh))
    UCRP["CAold"][OH]=np.multiply(UCoul["Aold"][OH],crp[OH])
    UC2["CAold"][OH]=np.multiply(UCoul["Aold"][OH],c2[OH])
    # New Acceptor
    if newaccept != -1:
        Ulj["Anew"].append(calc_lj(target,newaccept,eps,sig))
        UCoul["Anew"].append(calc_coul(target,newaccept,qo,qh))
    else:
        Ulj["Anew"].append(0.0)
        UCoul["Anew"].append(0.0)
    UCRP["CAnew"][OH]=np.multiply(UCoul["Anew"][OH],crp[OH])
    UCRP["LJAnew"][OH]=np.multiply(Ulj["Anew"][OH],crp[OH])
    UC2["CAnew"][OH]=np.multiply(UCoul["Anew"][OH],c2[OH])
    UC2["LJAnew"][OH]=np.multiply(Ulj["Anew"][OH],c2[OH])

    tmpljd=0.0
    # Donors
    Ulj["D"].append(0.0)
    UCoul["D"].append(0.0)
    for item in [donor1,donor2,donor3]:
        if item != -1:
            Ulj["D"][OH]   += calc_lj(target,item,eps,sig)
            UCoul["D"][OH] += calc_coul(target,donor1,qo,qh)
    UCRP["LJD"][OH]=np.multiply(Ulj["D"][OH],crp[OH])
    UCRP["CD"][OH]=np.multiply(UCoul["D"][OH],crp[OH])
    UC2["LJD"][OH]=np.multiply(Ulj["D"][OH],c2[OH])
    UC2["CD"][OH]=np.multiply(UCoul["D"][OH],c2[OH])

avcrp,avc2 = {},{}
for key in UCRP:
    avcrp[key]=np.sum(UCRP[key],axis=0)
    np.savetxt(key+'dcrp_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps,avcrp[key]],fmt="%2.5f")

for key in UC2:
    avc2[key]=np.divide(np.sum(UC2[key],axis=0),norm_t,out=np.zeros_like(np.sum(UC2[key],axis=0)),where=norm_t!=0)
    np.savetxt(key+'dframe_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps,avc2[key]],fmt="%2.5f")
"""

# This calulates the information for the jump angle and histograms it
tht,bedge = np.histogram(theta, bins=100,range=(0,3))
split=3/100.0
bins=[]
print(len(tht),len(bins))
for b in range(100):
    bins.append(b*split+split/2)

# Defines the step size
steps = np.linspace(0,ntimes-1,num=ntimes)
# Defines the normalization at time 0
norm_t[0]=len(OHs)

# Sets the special c1s
for i in range(len(norm_t)):
    n_hb = norm_t[i]
    if n_hb >= 5:
        c1_1[i],c1_2[i],c1_3[i],c1_4[i],c1_5[i] = C1[i],C1[i],C1[i],C1[i],C1[i]
    elif n_hb == 4:
        c1_1[i],c1_2[i],c1_3[i],c1_4[i],c1_5[i] = C1[i],C1[i],C1[i],C1[i], 0.0
    elif n_hb == 3:
        c1_1[i],c1_2[i],c1_3[i],c1_4[i],c1_5[i] = C1[i],C1[i],C1[i],0.0, 0.0
    elif n_hb == 2:
        c1_1[i],c1_2[i],c1_3[i],c1_4[i],c1_5[i] = C1[i],C1[i],0.0,0.0,0.0
    elif n_hb == 1:
        c1_1[i],c1_2[i],c1_3[i],c1_4[i],c1_5[i] = C1[i], 0.0, 0.0, 0.0, 0.0
    else: 
        c1_1[i],c1_2[i],c1_3[i],c1_4[i],c1_5[i] = 0.0, 0.0, 0.0, 0.0, 0.0
    


# Outputs the oo distance
np.savetxt('oo_dist.dat', np.c_[oo_dist])

# Outputs all the calculated correlation functions
#np.savetxt('ec1_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps,EC1], fmt="%2.5f")
#np.savetxt('lc1_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps,LC1], fmt="%2.5f")
#np.savetxt('theta_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[bins,tht], fmt="%2.5f")
#np.savetxt('crp_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, CRP], fmt="%2.5f")
#np.savetxt('hbnd_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, norm_t/norm_t[0]], fmt="%2.5f")
#np.savetxt('framec1_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, C1, norm_t], fmt="%2.5f")
#np.savetxt('framec2_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, C2, norm_t], fmt="%2.5f")
#np.savetxt('framec3_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, C3, norm_t], fmt="%2.5f")
#np.savetxt('framech1_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, CH1, norm_t], fmt="%2.5f")
#np.savetxt('framech2_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, CH2, norm_t], fmt="%2.5f")
#np.savetxt('framech3_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, CH3, norm_t], fmt="%2.5f")

# This part pulls the energies from the log file
print("Calculating the new energies")
edata = parse_log("log.lammps")
ke=edata["KinEng"][jindex]
lj=edata["E_vdwl"][jindex]
elec=np.add(edata["E_coul"][jindex],edata["E_long"][jindex])
tot=edata["TotEng"][jindex]
print("Finished the new energies")



import pickle


pickle.dump(jindex,open('jindex_'+str(nfile)+'_'+str(mol_name)+'.pckl','wb'))
pickle.dump(ke,open('jke_'+str(nfile)+'_'+str(mol_name)+'.pckl','wb'))
pickle.dump(lj,open('jlj_'+str(nfile)+'_'+str(mol_name)+'.pckl','wb'))
pickle.dump(elec,open('jelec_'+str(nfile)+'_'+str(mol_name)+'.pckl','wb'))
pickle.dump(tot,open('jtot_'+str(nfile)+'_'+str(mol_name)+'.pckl','wb'))
pickle.dump(theta,open('theta_'+str(nfile)+'_'+str(mol_name)+'.pckl','wb'))
pickle.dump(CRP,open('crp_'+str(nfile)+'_'+str(mol_name)+'.pckl', 'wb'))
pickle.dump(C1,open('framec1_'+str(nfile)+'_'+str(mol_name)+'.pckl', 'wb'))
pickle.dump(C2,open('framec2_'+str(nfile)+'_'+str(mol_name)+'.pckl', 'wb'))
pickle.dump(C3,open('framec3_'+str(nfile)+'_'+str(mol_name)+'.pckl', 'wb'))
pickle.dump(norm_t,open('norm_'+str(nfile)+'_'+str(mol_name)+'.pckl', 'wb'))

"""
np.savetxt('LJAold_init.out', np.c_[np.sum(Ulj["Aold"])])
np.savetxt('LJAnew_init.out', np.c_[np.sum(Ulj["Anew"])])
np.savetxt('LJD_init.out', np.c_[np.sum(Ulj["D"])])
np.savetxt('CAold_init.out', np.c_[np.sum(UCoul["Aold"])])
np.savetxt('CAnew_init.out', np.c_[np.sum(UCoul["Anew"])])
np.savetxt('CD_init.out', np.c_[np.sum(UCoul["D"])])
"""

print("End time %2.5f" % (time.time()-tread))                
