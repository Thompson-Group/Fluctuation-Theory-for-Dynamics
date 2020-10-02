#!/usr/bin/env python
# This code is used to calculate the activation energies of diffusion and reorientation based on the fluctuations method.
# It should run for either the NVT or NPT ensembles.

from scipy import stats
from scipy.optimize import curve_fit
from scipy.integrate import quad
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from argparse import RawTextHelpFormatter

def Corr_Fit(t,m, b):
    return m*t+b

def Error(value):
    err_value = np.array(value).std(0) * t_val
    return err_value

def E_Mult_Prop(a,ea,b,eb):
    return a*b*np.sqrt((ea/a)**2+(eb/b)**2)

def Order1_Prediction(x,T,M,dM):
    bzero=(1.0/(kb*T))
    return (1.0/6.0)*(M + dM*(x-bzero))*conv_d

def Order2_Prediction(x,T,M,dM,d2M):
    bzero=(1.0/(kb*T))
    return (1.0/6.0)*(M + dM*(x-bzero) + 0.5*d2M*(x-bzero)**2)*conv_d

def Order3_Prediction(x,T,M,dM,d2M,d3M):
    bzero=(1.0/(kb*T))
    return (1.0/6.0)*(M+dM*(x-bzero) + 0.5*d2M*(x-bzero)**2+(1.0/6.0)*d3M*(x-bzero)**3)*conv_d

def Order4_Prediction(x,T,M,dM,d2M,d3M,d4M):
    bzero=(1.0/(kb*T))
    return (1.0/6.0)*(M - dM*(x-bzero) + 0.5*d2M*(x-bzero)**2+(1.0/6.0)*d3M*(x-bzero)**3+(1.0/24.0)*d4M*(x-bzero)**4)*conv_d

# Read in Arguments
parser = argparse.ArgumentParser(description='''Calculates the msd fit''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
parser.add_argument('-blocks', help="Number of blocks")
parser.add_argument('-mol', help ="Molecule Name")
parser.add_argument('-temp', help = "Temperature (in Kelvin)")
parser.add_argument('-skip', help = "Data to skip")
args = parser.parse_args()

inputfile = str(args.inp)
nblocks = int(args.blocks)
mol_name = str(args.mol)
T=float(args.temp)
skip=int(args.skip)

#Physical Constants
kb=0.0019872041
conv_d=(1.0E12)/(1.0E16)*(1.0E5)
# Read in input file
inp_names=np.genfromtxt(inputfile, usecols=0, dtype=(str),unpack=True)
print(inp_names)

# Calculate important quantities
t_val=stats.t.ppf(0.95,nblocks-1)/np.sqrt(nblocks)


# Block Initializations
m_bl=[[] for x in range(0,nblocks)]
dm_bl=[[] for x in range(0,nblocks)]
d2m_bl=[[] for x in range(0,nblocks)]
d3m_bl=[[] for x in range(0,nblocks)]
d4m_bl=[[] for x in range(0,nblocks)]
b_bl=[[] for x in range(0,nblocks)]
db_bl=[[] for x in range(0,nblocks)]
d2b_bl=[[] for x in range(0,nblocks)]
d3b_bl=[[] for x in range(0,nblocks)]
d4b_bl=[[] for x in range(0,nblocks)]

beta=[[] for x in range(0,200)]
order1pred=[[] for x in range(0,200)]
order2pred=[[] for x in range(0,200)]
order3pred=[[] for x in range(0,200)]
order4pred=[[] for x in range(0,200)]
order1prederr=[[] for x in range(0,200)]
order2prederr=[[] for x in range(0,200)]
order3prederr=[[] for x in range(0,200)]
order4prederr=[[] for x in range(0,200)]
order1pred_bl=[[] for x in range(0,nblocks)]
order2pred_bl=[[] for x in range(0,nblocks)]
order3pred_bl=[[] for x in range(0,nblocks)]
order4pred_bl=[[] for x in range(0,nblocks)]

# Initialize Block Array

iindex=0
jindex=0

for item1 in inp_names:
    jindex=0
    for item2 in inp_names:
        if jindex >= iindex:
            item3=item2
            item4=item2
            for block in range(0,nblocks):              
                time, cab = np.genfromtxt('bl_' + str(block) + '_' + str(mol_name) + '_msd.dat', usecols=(0,1), unpack=True)
                dcab = np.genfromtxt('bl_' + str(block) + '_' + item1 + '_' + str(mol_name) + '_msd.dat', usecols=(1), unpack=True)
                time=time[skip:]
                cab=cab[skip:]
                dcab=dcab[skip:]
                popt_cab, pcov_cab = curve_fit(Corr_Fit, time, cab, p0=(1,0))
                m_bl[block]=popt_cab[0]
                b_bl[block]=popt_cab[1]
                print("Block %s First Derivative:" % block)
                popt_dcab, pcov_dcab = curve_fit(Corr_Fit,time, dcab, p0=(1,0))
                dm_bl[block]=popt_dcab[0]
                db_bl[block]=popt_dcab[1]
                print("Block %s Second Derivative:" % block)
                d2cab =np.genfromtxt('bl_' + str(block) + '_' + item1 + '_' + item2 + '_' + str(mol_name) + '_msd.dat', usecols=(1), unpack=True)
                d2cab=d2cab[skip:]
                popt_d2cab, pcov_d2cab = curve_fit(Corr_Fit, time, d2cab, p0=(1,0))
                d2m_bl[block]=popt_d2cab[0]
                d2b_bl[block]=popt_d2cab[1]
                print("Block %s Third Derivative:" % block)
                d3cab = np.genfromtxt('bl_'+str(block)+'_'+item1+'_'+item2+'_'+item3+'_'+mol_name+'_msd.dat', usecols=(1), unpack=True)
                d3cab=d3cab[skip:]
                popt_d3cab, pcov_d3cab = curve_fit(Corr_Fit, time, d3cab, p0=(1,0))
                d3m_bl[block]=popt_d3cab[0]
                d3b_bl[block]=popt_d3cab[1]
                print("Block %s Fourth Derivative:" % block)
                d4cab = np.genfromtxt('bl_'+str(block)+'_'+item1+'_'+item2+'_'+item3+'_'+item4+'_'+mol_name+'_msd.dat', usecols=(1), unpack=True)
                d4cab=d4cab[skip:]
                popt_d4cab, pcov_d4cab = curve_fit(Corr_Fit, time, d4cab, p0=(1,0))
                d4m_bl[block]=popt_d4cab[0]
                d4b_bl[block]=popt_d4cab[1]
                

                for b in range(0,200):
                    beta[b]=0.5+float(b)*0.01
                    order1pred[b]=Order1_Prediction(beta[b],T,m_bl[block], dm_bl[block])
                    order2pred[b]=Order2_Prediction(beta[b],T,m_bl[block], dm_bl[block], d2m_bl[block])
                    order3pred[b]=Order3_Prediction(beta[b],T,m_bl[block], dm_bl[block], d2m_bl[block], d3m_bl[block])
                    order4pred[b]=Order4_Prediction(beta[b],T,m_bl[block], dm_bl[block], d2m_bl[block], d3m_bl[block], d4m_bl[block])
                order1pred_bl[block]=order1pred
                order2pred_bl[block]=order2pred
                order3pred_bl[block]=order3pred
                order4pred_bl[block]=order4pred
                
            # I apologize for the lazy coding in the next section - for clarity:
            # Yes I am intentionally writing all the errors out to files just so they
            # can be read in later in the same program during the total calculation 
            err_m = Error(m_bl)
            err_dm = Error(dm_bl)
            err_d2m = Error(d2m_bl)
            err_d3m = Error(d3m_bl)
            err_d4m = Error(d4m_bl)
            err_b = Error(b_bl)
            err_db = Error(db_bl)
            err_d2b = Error(d2b_bl)
            err_d3b = Error(d3b_bl)
            err_d4b = Error(d4b_bl)
            
            order1prederr=np.array(order1pred_bl).std(0)
            order2prederr=np.array(order2pred_bl).std(0)
            order3prederr=np.array(order3pred_bl).std(0)
            order4prederr=np.array(order4pred_bl).std(0)
            order1prederr=[x * t_val for x in order1prederr]
            order2prederr=[x * t_val for x in order2prederr]
            order3prederr=[x * t_val for x in order3prederr]
            order4prederr=[x * t_val for x in order4prederr]

            np.savetxt("fiterr_"+item1+'_'+item2+'_'+item3+'_'+item4+'_'+mol_name+'_msd.dat', np.c_[err_m,err_b, err_dm,err_db, err_d2m,err_d2b, err_d3m,err_d3b, err_d4m, err_d4b], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+item4+'_'+str(mol_name)+'_errmsd.dat', np.c_[order1prederr, order2prederr, order3prederr, order4prederr], fmt='%s')


        jindex+=1
    iindex+=1

print("Beginning total calculation")
# Begin total Calculation
iindex=0
jindex=0
for item1 in inp_names:
    jindex=0
    for item2 in inp_names:
        if jindex >= iindex:
            item3=item2
            item4=item2
            time, cab = np.genfromtxt(str(mol_name) + '_msd.dat', usecols=(0,1,), unpack=True)
            dcab = np.genfromtxt(item1 + '_' + str(mol_name) + '_msd.dat', usecols=(1), unpack=True)
            time=time[skip:]
            cab=cab[skip:]
            dcab=dcab[skip:]
            popt_cab, pcov_cab = curve_fit(Corr_Fit, time, cab, p0=(.2,.3))
            m=popt_cab[0]
            b=popt_cab[1]
            print("First Total Derivative:")
            popt_dcab, pcov_dcab = curve_fit(Corr_Fit,time, dcab, p0=(0.2,0.3))
            dm=popt_dcab[0]
            db=popt_dcab[1]
            print("Second Total Derivative:")
            d2cab =np.genfromtxt(item1 + '_' + item2 + '_' + str(mol_name) + '_msd.dat', usecols=(1), unpack=True)
            d2cab=d2cab[skip:]
            popt_d2cab, pcov_d2cab = curve_fit(Corr_Fit, time, d2cab, p0=(0.2,0.3))
            d2m=popt_d2cab[0]
            d2b=popt_d2cab[1]
            print("Third Total Derivative:")
            d3cab = np.genfromtxt(item1+'_'+item2+'_'+item3+'_'+mol_name+'_msd.dat', usecols=(1), unpack=True)
            d3cab=d3cab[skip:]
            popt_d3cab, pcov_d3cab = curve_fit(Corr_Fit, time, d3cab, p0=(0.2,0.3))
            d3m=popt_d3cab[0]
            d3b=popt_d3cab[1]
            print("Fourth Total Derivative")
            d4cab = np.genfromtxt(item1+'_'+item2+'_'+item3+'_'+item4+'_'+mol_name+'_msd.dat', usecols=(1), unpack=True)
            d4cab=d4cab[skip:]
            popt_d4cab, pcov_d4cab = curve_fit(Corr_Fit, time, d4cab, p0=(1,0))
            d4m=popt_d4cab[0]
            d4b=popt_d4cab[1]

            err_m, err_b, err_dm, err_db, err_d2m, err_d2b, err_d3m, err_d3b, err_d4m, err_d4b= np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+item4+'_'+str(mol_name)+'_msd.dat', usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
            order1prederr,order2prederr,order3prederr,order4prederr = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+item4+'_'+str(mol_name)+'_errmsd.dat', usecols=(0,1,2,3), unpack=True)
            for bi in range(0,200):
                beta[bi]=0.5+float(bi)*0.01
                order1pred[bi]=Order1_Prediction(beta[bi],T,m,dm)
                order2pred[bi]=Order2_Prediction(beta[bi],T,m, dm, d2m)
                order3pred[bi]=Order3_Prediction(beta[bi],T,m, dm, d2m, d3m)
                order4pred[bi]=Order4_Prediction(beta[bi],T,m, dm, d2m, d3m, d4m)
                

            filepath=str('fit_results_'+item1+'_'+item2+'_'+str(mol_name)+'_msd.dat')
            fout=open(filepath,'w')
            fout.write("Correlation Function: MSD(t)\n")
            fout.write("Item 1: %s Item 2: %s\n" % (item1, item2))
            fout.write("  m: %s   err: %s\n" % (m,err_m))
            fout.write("  b: %s   err: %s\n" % (b, err_b))
            fout.write("  dm: %s   err: %s\n" % (dm, err_dm))
            fout.write("  db: %s   err: %s\n" % (db, err_db))
            fout.write("  d2m: %s   err: %s\n" % (d2m, err_d2m))
            fout.write("  d2b: %s   err: %s\n" % (d2b, err_d2b))
            fout.write(" d3m: %s   err: %s\n" % (d3m, err_d3m))
            fout.write(" d3b: %s   err: %s\n" % (d3b, err_d3b))
            fout.write(" d4m: %s   err: %s\n" % (d4m, err_d4m))
            fout.write(" d4b: %s   err: %s\n" % (d4b, err_d4b))
            ea = -dm/m
            err_ea = (dm/m)*np.sqrt((err_dm/dm)**2+(err_m/m)**2)
            fout.write(" Ea: %s    err: %s\n" % (ea, err_ea))
            fout.close()
            cabfit=Corr_Fit(time,*popt_cab)
            dcabfit=Corr_Fit(time, *popt_dcab)
            d2cabfit=Corr_Fit(time, *popt_d2cab)
            d3cabfit=Corr_Fit(time, *popt_d3cab)
            d4cabfit=Corr_Fit(time, *popt_d4cab)
            print(len(order2prederr))
            np.savetxt('total_fit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_msd.dat', np.c_[time,cabfit, dcabfit,d2cabfit], fmt='%s')
            np.savetxt('total_fit_'+str(item1)+'_'+str(item2)+'_'+str(item3)+'_'+str(mol_name)+'_msd.dat', np.c_[time, d3cabfit], fmt='%s')
            np.savetxt('total_fit_'+str(item1)+'_'+str(item2)+'_'+str(item3)+'_'+str(item4)+'_'+str(mol_name)+'_msd.dat', np.c_[time, d4cabfit], fmt='%s')
            np.savetxt('int_fit_'+str(item1)+'_'+str(mol_name)+'_msd.dat', np.c_[beta, order1pred, order1prederr], fmt='%s')
            np.savetxt('int_fit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_msd.dat', np.c_[beta, order2pred, order2prederr], fmt='%s')
            np.savetxt('int_fit_'+str(item1)+'_'+str(item2)+'_'+str(item3)+'_'+str(mol_name)+'_msd.dat', np.c_[beta, order3pred, order3prederr], fmt='%s')
            np.savetxt('int_fit_'+str(item1)+'_'+str(item2)+'_'+str(item3)+'_'+item4+'_'+str(mol_name)+'_msd.dat', np.c_[beta, order4pred, order4prederr], fmt='%s')
        jindex+=1
    iindex+=1
