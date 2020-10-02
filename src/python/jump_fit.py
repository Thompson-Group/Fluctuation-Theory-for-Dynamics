#!/usr/bin/env python
import numpy as np
from scipy.optimize import curve_fit
from numpy import exp, sin, sqrt
from scipy import stats
import sys,os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-corr', default='crp',type=str,help="Correlation function to fit, can be crp or frame")
parser.add_argument('-blocks', default=5, type=int, help="Number of blocks - note must match actual number!")
parser.add_argument('-cstart', default=0, type=int, help="Number of points to skip at start of fit")
parser.add_argument('-cend', default=5000, type=int, help="Number of points to skip at the end of fit")
parser.add_argument('-dcstart', default=0, type=int, help="Number of points to skip at start of fit")
parser.add_argument('-dcend', default=5000, type=int, help="Number of points to skip at the end of fit")
parser.add_argument('-degree', default=2, type=int, help="Degree of the Legendre Polynomial (default = 2)")
parser.add_argument('-keyoverride', default='all', type=str, help="Set if you want only one component")
parser.add_argument('-tidy', default=0, type=int, help="Set to 0 to skip, set to 1 to do post calculations")
parser.add_argument('-T', default=298.15, type=float, help="Temperature of the run")
parser.add_argument('-supress_output', default=0, type=int, help="Set to 1 to show no output, 0 shows all output")
args = parser.parse_args()

corr_func = args.corr
nblocks   = args.blocks
cstart    = args.cstart
cend      = args.cend
dcstart    = args.dcstart
dcend      = args.dcend
nl        = args.degree
setkey    = args.keyoverride
tidy      = args.tidy
out       = args.supress_output
T         = args.T

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

def single_exponential(x,A1, k1):
    return A1*exp(-k1*x)

def single_deriv_exponential(x, A, k, dA1, dk1):
    return (dA1-x*dk1*A[0])*exp(-k[0]*x)

def double_exponential(x, A1, k1, k2):
    return A1*exp(-k1*x)+(1-A1)*exp(-k2*x)

def double_deriv_exponential(x, A, k, dA1, dk1, dk2):
    return (dA1-x*dk1*A[0])*exp(-k[0]*x)+(-dA1-x*dk2*A[1])*exp(-k[1]*x)

def triple_exponential(x, A1,A2, k1, k2,k3):
    return A1*exp(-k1*x)+A2*exp(-k2*x)+(1-A1-A2)*exp(-k3*x)

def triple_deriv_exponential(x, A, k, dA1,dA2, dk1, dk2,dk3):
    return (dA1-x*dk1*A[0])*exp(-k[0]*x)+(dA2-x*dk2*A[1])*exp(-k[1]*x)+(-dA1-dA2-x*dk3*A[2])*exp(-k[2]*x)

def single_parse_exp_popt(popt):
    A,k = [],[]
    A=[popt[0]]
    k=[popt[1]]
    return A,k

def single_parse_dexp_popt(popt):
    dA,dk = [popt[0]],[popt[1]]
    return dA, dk

def double_parse_exp_popt(popt):
    A,k = [],[]
    if popt[1] > popt[2]:
        A=[popt[0],1-popt[0]]
        k=[popt[1],popt[2]]
    else:
        A=[1-popt[0],popt[0]]
        k=[popt[2],popt[1]]
    return A,k

def double_parse_dexp_popt(popt):
    dA,dk = [popt[0],-popt[0]],[popt[1],popt[2]]
    return dA, dk

def triple_parse_exp_popt(popt):
    A,k = [],[] # A is first 0,1 k is 2,3,4
    if popt[3] > popt[4]:
        if popt[2] > popt[3]: # 2 3 4
            A=[popt[0],popt[1],1-popt[0]-popt[1]]
            k=[popt[2],popt[3],popt[4]]
        elif popt[2] < popt[4]: # 3 4 2
            A=[popt[1],1-popt[0]-popt[1],popt[0]]
            k=[popt[3],popt[4],popt[2]]
        else: # 3 2 4
            A=[popt[1],popt[0],1-popt[0]-popt[1]]
            k=[popt[3],popt[2],popt[4]]
    else:
        if popt[2] < popt[3]: # 4 3 2 
            A=[1-popt[0]-popt[1],popt[1],popt[0]]
            k=[popt[4],popt[3],popt[2]]
        elif popt[2] > popt[4]: # 2 4 3
            A=[popt[0],1-popt[0]-popt[1],popt[1]]
            k=[popt[2],popt[4],popt[3]]
        else: # 4 2 3
            A=[1-popt[0]-popt[1], popt[0],popt[1]]
            k=[popt[4],popt[2],popt[3]]
    return A,k

def triple_parse_dexp_popt(popt):
    dA,dk = [popt[0],popt[1],-popt[0]-popt[1]],[popt[2],popt[3],popt[4]]
    return dA, dk

def of_of_theta(x,nl):
    #val = (1-sin(5/float(nl))/(5*sin(x/float(nl))))
    val = 1-(1/(2*nl+1))*sin((n+0.5)*x)/sin(x/2)
    return val

def do_cn_fit(xval, data, edata,bl_data, bl_edata):
    popt,pcov=curve_fit(triple_exponential,xval[cstart:cend],data[cstart:cend], p0=te,xtol=1e-6,ftol=1e-6)
    A,k = triple_parse_exp_popt(popt)
    bl_A, bl_k = [],[]
    bl_fit=[]
    for b in range(nblocks):
        popt,pcov=curve_fit(triple_exponential,xval[cstart:cend],bl_data[b][cstart:cend], p0=te)
        tmpA,tmpk = triple_parse_exp_popt(popt)
        bl_fit.append(triple_exponential(xval,tmpA[0],tmpA[1],tmpk[0],tmpk[1],tmpk[2]))
        bl_A.append(tmpA)
        bl_k.append(tmpk)
    np.savetxt(corr_func+'fit_results.dat', np.c_[xval,triple_exponential(xval, A[0],A[1], k[0],k[1],k[2]),np.std(bl_fit,axis=0)*t_val])
    err={}
    err["k"]=np.std(bl_k,axis=0)*t_val
    err["A"]=np.std(bl_A,axis=0)*t_val
    for key in edata:
        dpopt,dpcov=curve_fit(lambda x, dA1,dA2, dk1, dk2,dk3: triple_deriv_exponential(x, A, k, dA1,dA2, dk1, dk2,dk3), xval[dcstart:dcend], edata[key][dcstart:dcend],p0=dte,xtol=1e-6,ftol=1e-6)
        dA, dk = triple_parse_dexp_popt(dpopt)
        bl_dA,bl_dk,bl_ea = [],[],[]
        bl_efit = []
        for b in range(nblocks):
            dpopt,dpcov=curve_fit(lambda x, dA1,dA2, dk1, dk2,dk3: triple_deriv_exponential(x, A, k, dA1,dA2, dk1, dk2,dk3), xval[dcstart:dcend], bl_edata[key][b][dcstart:dcend],p0=dte,xtol=1e-6,ftol=1e-6)
            tmpdA,tmpdk = triple_parse_dexp_popt(dpopt)
            bl_efit.append(triple_deriv_exponential(xval,A,k,tmpdA[0],tmpdA[1],tmpdk[0],tmpdk[1],tmpdk[2]))
            bl_dA.append(tmpdA)
            bl_dk.append(tmpdk)
        bl_ea = -np.divide(bl_dk,bl_k)
        err["dk"]=np.std(bl_dk,axis=0)*t_val
        err["dA"]=np.std(bl_dA,axis=0)*t_val
        err["ea"]=np.std(bl_ea,axis=0)*t_val
        print_data(key,3,A,k,dA,dk,err)
        np.savetxt(corr_func+'fit_results_'+key+'.dat', np.c_[xval,triple_deriv_exponential(xval, A, k, dA[0],dA[1], dk[0], dk[1],dk[2]),np.std(bl_efit,axis=0)*t_val])
    return

def do_crp_fit(xval, data, edata,bl_data, bl_edata):
    popt,pcov=curve_fit(double_exponential,xval[cstart:cend],1-data[cstart:cend], p0=de,xtol=1e-6,ftol=1e-6)
    A,k = double_parse_exp_popt(popt)
    bl_A, bl_k = [],[]
    bl_fit=[]
    for b in range(nblocks):
        popt,pcov=curve_fit(double_exponential,xval[cstart:cend],1-bl_data[b][cstart:cend], p0=de)
        tmpA,tmpk = double_parse_exp_popt(popt)
        bl_fit.append(double_exponential(xval,tmpA[0],tmpk[0],tmpk[1]))
        bl_A.append(tmpA)
        bl_k.append(tmpk)
    np.savetxt('crpfit_results.dat', np.c_[xval,double_exponential(xval, A[0], k[0],k[1]),np.std(bl_fit,axis=0)*t_val])
    err={}
    err["k"]=np.std(bl_k,axis=0)*t_val
    err["A"]=np.std(bl_A,axis=0)*t_val 
    for key in edata:
        dpopt,dpcov=curve_fit(lambda x, dA1, dk1, dk2: double_deriv_exponential(x, A, k, dA1, dk1, dk2), xval[dcstart:dcend], -edata[key][dcstart:dcend],p0=dde,xtol=1e-6,ftol=1e-6)
        dA, dk = double_parse_dexp_popt(dpopt)
        bl_dA,bl_dk,bl_ea = [],[],[]
        bl_efit = []
        for b in range(nblocks):
            dpopt,dpcov=curve_fit(lambda x, dA1, dk1, dk2: double_deriv_exponential(x, A, k, dA1, dk1, dk2), xval[dcstart:dcend], -bl_edata[key][b][dcstart:dcend],p0=dde,xtol=1e-6,ftol=1e-6)
            tmpdA,tmpdk = double_parse_dexp_popt(dpopt)
            bl_efit.append(double_deriv_exponential(xval,A,k,tmpdA[0],tmpdk[0],tmpdk[1]))
            bl_dA.append(tmpdA)
            bl_dk.append(tmpdk)
        bl_ea = -np.divide(bl_dk,bl_k)
        err["dk"]=np.std(bl_dk,axis=0)*t_val
        err["dA"]=np.std(bl_dA,axis=0)*t_val
        err["ea"]=np.std(bl_ea,axis=0)*t_val
        print_data(key,2,A,k,dA,dk,err)
        np.savetxt('crpfit_results_'+key+'.dat', np.c_[xval,double_deriv_exponential(xval, A, k, dA[0], dk[0], dk[1]),np.std(bl_efit,axis=0)*t_val])
    return

def do_frame_fit(xval, data, edata,bl_data, bl_edata):
    popt,pcov=curve_fit(double_exponential,xval[cstart:cend],data[cstart:cend], p0=de,xtol=1e-6,ftol=1e-6)
    A,k = double_parse_exp_popt(popt)
    bl_A, bl_k = [],[]
    bl_fit=[]
    for b in range(nblocks):
        popt,pcov=curve_fit(double_exponential,xval[cstart:cend],bl_data[b][cstart:cend], p0=de)
        tmpA,tmpk = double_parse_exp_popt(popt)
        bl_fit.append(double_exponential(xval,tmpA[0],tmpk[0],tmpk[1]))
        bl_A.append(tmpA)
        bl_k.append(tmpk)
    np.savetxt(corr_func+'fit_results.dat', np.c_[xval,double_exponential(xval, A[0], k[0],k[1]),np.std(bl_fit,axis=0)*t_val])
    err={}
    print(np.divide(1,bl_k))
    err["k"]=np.std(bl_k,axis=0)*t_val
    err["A"]=np.std(bl_A,axis=0)*t_val
    for key in edata:
        dpopt,dpcov=curve_fit(lambda x, dA1, dk1,dk2: double_deriv_exponential(x, A, k, dA1, dk1, dk2), xval[dcstart:dcend], edata[key][dcstart:dcend],p0=dde,xtol=1e-6,ftol=1e-6)
        dA, dk = double_parse_dexp_popt(dpopt)
        bl_dA,bl_dk,bl_ea = [],[],[]
        bl_efit = []
        for b in range(nblocks):
            dpopt,dpcov=curve_fit(lambda x, dA1, dk1, dk2: double_deriv_exponential(x, A, k, dA1, dk1,dk2), xval[dcstart:dcend], bl_edata[key][b][dcstart:dcend],p0=dde,xtol=1e-6,ftol=1e-6)
            tmpdA,tmpdk = double_parse_dexp_popt(dpopt)
            bl_efit.append(double_deriv_exponential(xval,A,k,tmpdA[0],tmpdk[0],tmpdk[1]))
            bl_dA.append(tmpdA)
            bl_dk.append(tmpdk)
        bl_ea = -np.divide(bl_dk,bl_k)
        err["dk"]=np.std(bl_dk,axis=0)*t_val
        err["dA"]=np.std(bl_dA,axis=0)*t_val
        err["ea"]=np.std(bl_ea,axis=0)*t_val
        print_data(key,2,A,k,dA,dk,err)
        np.savetxt(corr_func+'fit_results_'+key+'.dat', np.c_[xval,double_deriv_exponential(xval, A, k, dA[0], dk[0],dk[1]),np.std(bl_efit,axis=0)*t_val])
    return

def f_of_theta(x,nl):
    val = (1-sin((nl+0.5)*x)/((2*nl+1)*sin(x/2)))
    return val


def norm_theta(xval, theta):
    int_theta = np.trapz(theta,xval)
    norm_theta = np.divide(theta, int_theta)
    return norm_theta, int_theta

def do_theta_int(xval, data, edata,bl_data, bl_edata):
    print("nl is %s" % nl)
    theta, int_theta = norm_theta(xval, data)
    print("norm_const",int_theta)
    print("average angle: %s " % np.trapz(np.multiply(theta,xval)*57.2958, xval))

    avtheta = np.trapz(np.multiply(f_of_theta(xval,nl),theta),xval)
    bl_theta = np.divide(bl_data,int_theta)
    np.savetxt(corr_func+str(nl)+'_int_results.dat', np.c_[xval,theta,np.std(bl_theta,axis=0)*t_val])
    bl_avtheta = []
    bltmp = []
    for b in range(nblocks):
        bl_avtheta.append(np.trapz(np.multiply(f_of_theta(xval,nl),bl_theta[b]),xval))
        bltmp.append(np.trapz(np.multiply(bl_theta[b],xval)*57.2958, xval))
    err={}
    tmperr = np.std(bltmp,axis=0)*t_val
    print("average angle: %s +/- %s" % (np.trapz(np.multiply(theta,xval)*57.2958, xval),tmperr))
    err["avtheta"]=np.std(bl_avtheta,axis=0)*t_val
    for key in edata:
        avetheta = np.trapz(np.multiply(f_of_theta(xval,nl),edata[key]/int_theta),xval)
        bl_etheta = np.divide(bl_edata[key],int_theta)
        bl_avetheta = []
        for b in range(nblocks):
            bl_avetheta.append(np.trapz(np.multiply(f_of_theta(xval,nl),bl_edata[key][b]/int_theta),xval))
        bl_ea = np.divide(bl_avetheta,bl_avtheta,out=np.zeros_like(bl_avetheta),where=bl_avtheta!=0)
        err["avetheta"]=np.std(bl_avetheta,axis=0)*t_val
        err["ea"]=np.std(bl_ea,axis=0)*t_val
        print_data(key,0,0, avtheta,0, avetheta,err)
        np.savetxt(corr_func+str(nl)+'_int_results_'+key+'.dat', np.c_[xval,edata[key]/int_theta,np.std(bl_etheta,axis=0)*t_val])
    return



def print_data(item,n, A, k, dA, dk, err):
    if out == 0: print(color.BOLD + color.PURPLE +  color.UNDERLINE + "#TYPE KEY IS %s" % item + color.END)
    if out == 0:print(color.BOLD + "#Type  Unit  num: value err"+color.END)
    if corr_func != "theta" and out == 0:
        for i in range(n):
            print("t  (ps)          %d: % 09.5f % 09.5f" % (i+1,1/k[i], 1/k[i]*err["k"][i]/k[i]))
            print("k  (ps^-1)       %d: % 09.5f % 09.5f" % (i+1,k[i],err["k"][i]))
            print("dk (kcal/mol/ps) %d: % 09.5f % 09.5f" % (i+1,dk[i],err["dk"][i]))
            print("A  (Unitless)    %d: % 09.5f % 09.5f" % (i+1,A[i],err["A"][i]))
            print("dA (Unitless)    %d: % 09.5f % 09.5f" % (i+1,dA[i],err["dA"][i]))
            print(color.RED + color.BOLD+"Ea (kcal/mol)    %d: % 09.5f % 09.5f" % (i+1,-dk[i]/k[i],err["ea"][i]) + color.END)
            if item == "e":
                f=open('ea_%s_%d.dat' % (corr_func,i+1), 'w')
                f.write("Ea (kcal/mol)    %8s: % 09.10f % 09.10f TS % 09.10f % 09.10f\n" % (item,-dk[i]/k[i],err["ea"][i],k[i],err["k"][i]))
                f.close()
            else:
                f=open('ea_%s_%d.dat' % (corr_func,i+1), 'a')
                f.write("Ea (kcal/mol)    %8s: % 09.10f % 09.10f TS % 09.10f % 09.10f\n" % (item,-dk[i]/k[i],err["ea"][i],k[i],err["k"][i]))
                f.close()
    elif out == 0:
        print("thta  (deg)          %d: % 09.5f % 09.5f" % (1,k, err["avtheta"]))
        print("ethta (deg kcal/mol) %d: % 09.5f % 09.5f" % (1,dk, err["avetheta"]))
        print(color.RED + color.BOLD+"Ea (kcal/mol)        %d: % 09.5f % 09.5f" % (1,-dk/k,err["ea"]) + color.END)
        if item == "e":
            f=open('ea_%s%s.dat' % (corr_func,nl), 'w')
            f.write("Ea (kcal/mol)    %8s: % 09.10f % 09.10f TS % 09.10f % 09.10f\n" % (item,-dk/k,err["ea"],k,err["avtheta"]))
            f.close()
        else:
            f=open('ea_%s%s.dat' % (corr_func,nl), 'a')
            f.write("Ea (kcal/mol)    %8s: % 09.10f % 09.10f TS % 09.10f % 09.10f\n" % (item,-dk/k,err["ea"],k,err["avtheta"]))
            f.close()
    return

def propagate_div(q1,q2,eq1,eq2):
    err=np.divide(q1,q2,where=q2!=0, out=np.zeros_like(q1))*sqrt(np.divide(eq1,q1,where=q1!=0,out=np.zeros_like(eq1))**2.+np.divide(eq2,q2,where=q2!=0,out=np.zeros_like(eq2))**2.)
    return err


def propagate_mult(q1,q2,eq1,eq2):
    err=(q1*q2)*sqrt(np.divide(eq1,q1,where=q1!=0,out=np.zeros_like(eq1))**2.+np.divide(eq2,q2,where=q2!=0,out=np.zeros_like(eq2))**2.)
    return err

def propagate_add(eq1,eq2):
    err=sqrt(eq1**2.+eq2**2.)
    return err

def my_ceil(a, precision=0):
    return np.round(a + 0.5 * 10**(-precision), precision)

def my_floor(a, precision=0):
    return np.round(a - 0.5 * 10**(-precision), precision)

def zero_count(error):
    error = np.abs(error)
    return int(np.ceil(-np.log10(np.abs(error) - np.abs(np.floor(error)))))

def setprecision(value, error):
    prec = zero_count(error)
    if prec == 1:
        prec =2
    value = "%0.*f" % (prec,np.round(value,prec))
    error = int(my_ceil(error,prec)*10**prec)
    return value, np.abs(error)

def read_tidy():
    ea_theta,err_theta,ftheta,err_ftheta = np.zeros(len(inp_n)),np.zeros(len(inp_n)),np.zeros(len(inp_n)),np.zeros(len(inp_n))
    ea_c, err_c = np.zeros(len(inp_n)),np.zeros(len(inp_n))
    ea_frame, err_frame,kframe, err_kframe = np.zeros(len(inp_n)),np.zeros(len(inp_n)),np.zeros(len(inp_n)),np.zeros(len(inp_n))
    ea_ko, err_ea_ko,ko,err_ko = np.zeros(len(inp_n)),np.zeros(len(inp_n)),np.zeros(len(inp_n)),np.zeros(len(inp_n))
    k_c, err_k_c = np.zeros(len(inp_n)),np.zeros(len(inp_n))
    if os.path.isfile("ea_theta%d.dat"%nl):
        ea_theta,err_theta,ftheta,err_ftheta=np.genfromtxt("ea_theta%d.dat"%nl, usecols=(3,4,6,7),unpack=True)
    if os.path.isfile("ea_c%d_3.dat"%nl):
        ea_c, err_c=np.genfromtxt("ea_c%d_3.dat"%nl,usecols=(3,4),unpack=True)
    if os.path.isfile("ea_framec%d_2.dat"%nl):
        ea_frame, err_frame,kframe, err_kframe=np.genfromtxt("ea_framec%d_2.dat"%nl,usecols=(3,4,6,7),unpack=True)
    if os.path.isfile("ea_crp_2.dat"):
        ea_ko, err_ea_ko,ko,err_ko=np.genfromtxt("ea_crp_2.dat",usecols=(3,4,6,7),unpack=True)
    if os.path.isfile("ea_c%d_3.dat"%nl):
        k_c, err_k_c = np.genfromtxt("ea_c%d_3.dat"%nl,usecols=(6,7),unpack=True)
    ea_jump, err_ea_jump = ea_ko+ea_theta, propagate_add(err_ea_ko,err_theta)
    # Timescale Reorientation
    t_c, err_t_c = np.divide(1,k_c[0]),np.divide(1,k_c[0])*np.divide(err_k_c[0],k_c[0])
    # Angle Distribution
    ftheta,err_ftheta = ftheta[0],err_ftheta[0]
    # Frame 
    kframe, err_kframe = kframe[0], err_kframe[0]
    tframe, err_tframe = np.divide(1,kframe),np.divide(1,kframe)*np.divide(err_kframe,kframe)
    # Jump Time
    ko,err_ko = ko[0],err_ko[0]
    to, err_to = np.divide(1,ko),np.divide(1,ko)*np.divide(err_ko,ko)
    # jump cont
    kjump, err_kjump = np.multiply(ko,ftheta),propagate_mult(ko,ftheta,err_ko,err_ftheta)
    tjump, err_tjump = np.divide(1,kjump),np.divide(1,kjump)*np.divide(err_kjump,kjump)
    # EJM
    kbot, err_kbot = kjump + kframe, propagate_add(err_kjump,err_kframe)
    kjrat, err_kjrat = np.divide(kjump,kbot,where=kbot!=0,out=np.zeros_like(kjump)), propagate_div(kjump,kbot,err_kjump,err_kbot)
    kfrat, err_kfrat = np.divide(kframe,kbot,where=kbot!=0,out=np.zeros_like(kframe)), propagate_div(kframe,kbot,err_kframe,err_kbot)
    kejm, err_kejm = kbot, err_kbot
    tejm, err_tejm = np.divide(1,kejm),np.divide(1,kejm)*np.divide(err_kejm,kejm)
    ea_c2_pred, err_c2_pred = kjrat*ea_jump + kfrat*ea_frame, propagate_add(propagate_mult(kjrat,ea_jump,err_kjrat,err_ea_jump),propagate_mult(kfrat,ea_frame,err_kfrat,err_frame))

    print("TIDYING UP JUMP CALCULATIONS - Legendre Polynomial Degree: %d" %nl)
    print("fthta = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (ftheta,err_ftheta,0,0))
    if k_c[0] != 0: print("Kc    = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (k_c[0],err_k_c[0],t_c,err_t_c))
    if ko != 0: print("Ko    = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (ko,err_ko,to,ko))
    if kjump != 0: print("Kjump = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (kjump,err_kjump,tjump,tjump))
    if kframe != 0: print("Kfram = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (kframe,err_kframe,tframe,err_tframe))
    print("Kjrat = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (kjrat,err_kjrat,0,0))
    print("Kfrat = % 09.5f +/- % 09.5f (% 09.5f +/- % 09.5f)" % (kfrat,err_kfrat,0,0))
    print("PRINTING INDIVIDUAL ACTIVATION ENERGIES")
    print("%8s %9s (%9s) : %9s (%9s) : %9s (%9s) : %9s (%9s) : %9s (%9s)" % ("#ITEM", "C2ACT","Error","JUMP","Error","FRAME","Error","THETA","Error", "ko", "ERROR"))
    actvol_conv = 83.14462*T
    for i in range(len(ea_theta)):
        if inp_n[i] == "vol":
            ea_c[i],err_c[i] = ea_c[i] * actvol_conv,err_c[i] *actvol_conv
            ea_jump[i], err_ea_jump[i] = ea_jump[i] *actvol_conv, err_ea_jump[i]*actvol_conv
            ea_frame[i],err_frame[i] = ea_frame[i]* actvol_conv,err_frame[i]* actvol_conv
            ea_theta[i], err_theta[i] = ea_theta[i]* actvol_conv, err_theta[i]* actvol_conv
            ea_ko[i], err_ea_ko[i] = ea_ko[i]* actvol_conv, err_ea_ko[i]* actvol_conv
            ea_c2_pred[i], err_c2_pred[i] = ea_c2_pred[i] * actvol_conv, err_c2_pred[i] * actvol_conv
        print("%8s % 09.5f (% 09.5f) : % 09.5f (% 09.5f) : % 09.5f (% 09.5f) : % 09.5f (% 09.5f) : % 09.5f (% 09.5f)" % (inp_n[i], ea_c[i], err_c[i], ea_jump[i], err_ea_jump[i], ea_frame[i],err_frame[i], ea_theta[i], err_theta[i],ea_ko[i], err_ea_ko[i]))
    print("PRINTING Combined ACTIVATION ENERGIES")
    print("%8s %9s (%9s) = %9s (%9s) + %9s (%9s) = %9s (%9s)" % ("#ITEM", "C2ACT","Error","JUMP","Error","FRAME","Error","C2COMP","Error"))
    jumpcont,err_jumpcont=np.zeros(len(ea_theta)),np.zeros(len(ea_theta))
    framecont,err_framecont=np.zeros(len(ea_theta)),np.zeros(len(ea_theta))
    for i in range(len(ea_theta)):
        jumpcont[i],err_jumpcont[i] = kjrat*ea_jump[i],propagate_mult(kjrat,ea_jump,err_kjrat,err_ea_jump)[i]
        framecont[i], err_framecont[i] = kfrat*ea_frame[i],propagate_mult(kfrat,ea_frame,err_kfrat,err_frame)[i]
        print("%8s % 09.5f (% 09.5f) = % 09.5f (% 09.5f) + % 09.5f (% 09.5f) = % 09.5f (% 09.5f)" % (inp_n[i],ea_c[i],err_c[i],jumpcont[i],err_jumpcont[i],framecont[i],err_framecont[i],ea_c2_pred[i],err_c2_pred[i]))
    """
    # Redoes precision
    to,err_to = setprecision(to,err_to)
    ftheta,err_ftheta = setprecision(ftheta,err_ftheta)
    tjump,err_tjump = setprecision(tjump,err_tjump)
    tframe, err_tframe = setprecision(tframe,err_tframe)
    t_c,err_t_c = setprecision(t_c,err_t_c)
    tejm,err_tejm = setprecision(tejm,err_tejm)
    print("TIMESCALES TABLE - LATEX")
    print("Order &  $\\tau_0$ & $\\overline{w}_n$ &  $\\tau_n^{jump}$ & $\\tau_n^{frame}$ & $\\tau_n$ & $\\tau_{n}^{EJM}$  \\\\")
    print("$n=%d$ & %s$_{%d}$ & %s$_{%d}$ & %s$_{%d}$ & %s$_{%d}$ & %s$_{%d}$ & %s$_{%d}$ \\\\" % (nl,to,err_to,ftheta,err_ftheta,tjump,err_tjump,tframe,err_tframe,t_c,err_t_c,tejm,err_tejm))
    print("END TIMESCALES TABLE ENTRY")

    print("Activation Energy Table Latex")
    print("Component &  $E_{a,0}$ & $E_{a,n}^{\\Delta\\theta}$ & $E_{a,n}^{jump}$ & $E_{a,n}^{frame}$ \\\\")
    for i in range(len(ea_theta)):
        ea_c[i], err_c[i] = setprecision(ea_c[i],err_c[i])
        ea_jump[i],err_ea_jump[i] = setprecision(ea_jump[i],err_ea_jump[i])
        ea_frame[i],err_frame[i]  = setprecision(ea_frame[i],err_frame[i])
        ea_theta[i],err_theta[i] = setprecision(ea_theta[i],err_theta[i])
        ea_ko[i],err_ea_ko[i] = setprecision(ea_ko[i],err_ea_ko[i])
        print("%s & %s$_{%d}$ & %s$_{%d}$ & %s$_{%d}$ & %s$_{%d}$\\\\" %  (inp_n[i],ea_ko[i],err_ea_ko[i],ea_theta[i],err_theta[i],ea_jump[i],err_ea_jump[i],ea_frame[i],err_frame[i]))
    print("Finish Act E table")
    print("Combination table latex")
    print("Component & $E_{a,n}$ & = & Jump  & + & Frame   & = & $E_{a,n}^{EJM}$  \\\\ & Calculated & & Contribution & & Contribution & & Predicted \\\\")
    for i in range(len(ea_theta)):
        ea_c2_pred[i],err_c2_pred[i] = setprecision(ea_c2_pred[i],err_c2_pred[i])
        jumpcont[i],err_jumpcont[i] = setprecision(jumpcont[i],err_jumpcont[i])
        framecont[i],err_framecont[i] = setprecision(framecont[i],err_framecont[i])
        print("%s & %s$_{%d}$ & = & %s$_{%d}$ & + & %s$_{%d}$ & = & %s$_{%d}$ \\\\" % (inp_n[i],ea_c[i],err_c[i],jumpcont[i], err_jumpcont[i], framecont[i], err_framecont[i],ea_c2_pred[i],err_c2_pred[i]))
    print("finish combination table latex")
    """
    print("JUMP CALCULATIONS COMPLETE")
    return

def do_tidy():
    read_tidy()
    return
    
    
t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)


inp_n = np.genfromtxt("flucts.inp", usecols=(0),dtype=str,unpack=True)
if setkey != 'all':
        inp_n = setkey

if tidy == 1:
    if nl == 0:
        for nl in range(1,4):
            do_tidy()
    else:
        do_tidy()
    sys.exit("Tidying up completed")





#initial_guesses
se=(0.9,0.3)
dse=(0.7,4)
de=(0.9, 0.3, 54.0)
dde=(0.7, 30, -0.9)
te=(0.05,0.2, 40, 2,0.3)
dte=(0.02,0.1, 3, 1,0.4)

print("JUMP CALCULATIONS BEGIN")

if corr_func not in ["framec1","framec2","framec3","crp","theta","c1","c2","c3"]:
    sys.exit("Error: That corr_func (%s) is as defined as 5/0" % corr_func)
subdir="OUT/"
xval, data = np.genfromtxt(subdir+"water_%s.dat" % corr_func,usecols=(0,1),unpack=True)
bl_data = []
for b in range(nblocks):
    bl_data.append(np.genfromtxt(subdir+"bl_%d_water_%s.dat" % (b,corr_func),usecols=(1),unpack=True))


pderiv_conv=138.0649*T # kbT in bar A^3 

edata,bl_edata={},{}
for item in inp_n:
    edata[item]=np.genfromtxt(subdir+"%s_water_%s.dat" % (item,corr_func),usecols=(1),unpack=True)
    bl_edata[item]=[]
    for b in range(nblocks):
        bl_edata[item].append(np.genfromtxt(subdir+"bl_%d_%s_water_%s.dat" % (b,item,corr_func),usecols=(1),unpack=True))
    if item == "vol":
        edata[item] = np.divide(edata[item],pderiv_conv)
        bl_edata[item] = np.divide(bl_edata[item],pderiv_conv)
if corr_func == "crp":
    do_crp_fit(xval,data,edata,bl_data,bl_edata)
elif corr_func in ["c1","c2","c3"]:
    do_cn_fit(xval,data,edata,bl_data,bl_edata)
elif corr_func in ["framec1","framec2","framec3"]:
    do_frame_fit(xval,data,edata,bl_data,bl_edata)
elif corr_func == "theta":
    do_theta_int(xval,data,edata,bl_data,bl_edata)
else:
    sys.exit("Error: That corr_func (%s) is as defined as 5/0" % corr_func)

print("JUMP CALCULATIONS COMPLETE")
