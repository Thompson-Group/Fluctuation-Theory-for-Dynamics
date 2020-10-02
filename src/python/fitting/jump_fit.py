#!usr/bin/env python
import numpy as np
from scipy.optimize import curve_fit
from numpy import exp

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def double_exponential(x, A1, k1, k2):
    return A1*exp(-k1*x)+(1-A1)*exp(-k2*x)

def double_deriv_exponential(x, A, k, dA1, dk1, dk2):
    return (dA1-x*dk1*A[0])*exp(-k[0]*x)+(-dA1-x*dk2*A[1])*exp(-k[1]*x)

def parse_exp_popt(popt):
    A,k = [],[]
    if popt[1] > popt[2]:
        A=[popt[0],1-popt[0]]
        k=[popt[1],popt[2]]
    else:
        A=[1-popt[0],popt[0]]
        k=[popt[2],popt[1]]
    return A,k

def parse_dexp_popt(popt):
    dA,dk = [popt[0],1-popt[0]],[popt[1],popt[2]]
    return dA, dk

def do_fit(xval, data, edata):
    popt,pcov=curve_fit(double_exponential,xval,1-data, p0=de)
    A,k = parse_exp_popt(popt)
    dpopt,dpcov=curve_fit(lambda x, dA1, dk1, dk2: double_deriv_exponential(x, A, k, dA1, dk1, dk2), xval, -edata,p0=dde)
    dA, dk = parse_dexp_popt(dpopt)
    print_data(A,k,dA,dk)
    return

def print_data(A, k, dA, dk):
    print(color.BOLD + color.PURPLE + "#Type    num: value err"+color.END)
    for i in range(2):
        print("Timescale  %d: %5.4f %s" % (i+1,k[i],0))
        print("dTimescale %d: %5.4f %s" % (i+1,dk[i],0))
        print(color.RED + color.BOLD+"Activ_Ener %d: %5.4f %s" % (i+1,-dk[i]/k[i],0) + color.END)
        print("A          %d: %5.4f %s" % (i+1,A[i],0))
        print("dA         %d: %5.4f %s" % (i+1,dA[i],0))
    return
    

#initial_guesses
de=(0.9, 0.3, 54.0)
dde=(0.7, 30, -0.9)

# Section to read flucts.inp

xval, data = np.genfromtxt("/panfs/pfs.local/scratch/thompson/e924p726/new_jump/test1/with-theta/water_crp.dat",usecols=(0,1),unpack=True)
edata=np.genfromtxt("/panfs/pfs.local/scratch/thompson/e924p726/new_jump/test1/with-theta/e_water_crp.dat",usecols=(1),unpack=True)

do_fit(xval,data,edata)

