#!/usr/bin/env python
import numpy as np
import pickle
import sys,os
import argparse
from scipy import stats
from read_input import user_input

"""
This script has two main roles:
    1) Parse all the correlation functions into something quickly readable
    2) Calculate the correlation function + its derivative.

These two options can be selected with the option flag, and should be run in that order.

Note this script has the support for the first and second derivatives included. Third and fourth derivatives are output as empty files currently, but will be added in at a later date. If these derivatives are needed, then one should use the old programs for calculating these init_segments.py and combine_segments.py.

This script is MUCH faster than its predecessors, especially when considering the ability to read binary files.
It also takes up much less storage space for the same set of calculations.
"""

def calc_avcab(cab, no):
    """
    Calculates the average, with support for a changing normalization
    """
    num = np.sum(cab,axis=0)
    norm = np.sum(no,axis=0)
    return np.divide(num,norm, where=norm!=0, out=np.zeros_like(num))

def zero_order(cab, no):
    """
    Calculates the zeroth order derivative - aka the correlation function!
    Does block averaging, and the error calculation
    """
    d0 = calc_avcab(cab,no)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d0 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d0.append(calc_avcab(cab[bstart:bend],no[bstart:bend]))
        np.savetxt(subdir+"bl_"+str(block)+"_"+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d0[block]])
    err_d0 = Error(bl_d0)
    np.savetxt(subdir+mol_name+'_'+corr_func+'.dat',np.c_[time,d0,err_d0])
    return

def zero_order_dist(cab):
    """
    Calculates the zeroth order derivative - aka the correlation function!
    Does block averaging, and the error calculation
    """
    notmp = np.ones(np.shape(cab))
    d0 = calc_avcab(cab,notmp)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d0 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d0.append(calc_avcab(cab[bstart:bend],notmp[bstart:bend]))
    err_d0 = Error(bl_d0)
    print("The average of the histogram values are: %13.9f +/- %13.9f" %(d0, err_d0))
    return

def calc_deriv(cab, encab, no):
    """ 
    Calculates the first derivative
    f' = - < dH f >
    """
    avencab=calc_avcab(encab,no)
    dencab = np.subtract(encab, avencab)
    dcab = -np.multiply(dencab,cab)
    #av_dcab = calc_avcab(dcab,nosplit*inputparam.segsplit)
    av_dcab = calc_avcab(dcab,no)
    return av_dcab

def first_order(cab,encab,no):
    """
    This does the first derivative and the blocking!
    """
    d1 = calc_deriv(cab, encab, no)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d1 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d1.append(calc_deriv(cab[bstart:bend],encab[bstart:bend],no[bstart:bend]))
        np.savetxt(subdir+"bl_"+str(block)+"_"+item+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d1[block]])
    err_d1 = Error(bl_d1)
    np.savetxt(subdir+item+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d1,err_d1])
    return 

def first_order_dist(cab,encab):
    """
    This does the first derivative and the blocking!
    """
    notmp = np.ones(np.shape(cab))
    d1 = calc_deriv(cab, encab, notmp)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d1 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d1.append(calc_deriv(cab[bstart:bend],encab[bstart:bend],notmp[bstart:bend]))
    err_d1 = Error(bl_d1)
    print("The derivative of the distribution is: %13.9f +/- %13.9f" % (d1,err_d1))
    return

def calc_deriv2(cab,encab1,encab2,no):
    """
    This calculates the second derivative as:
    f'' = <[dH^2 - <dH^2>]f>
    """
    #Second Derivative
    # <H1>, <H2>
    avencab1, avencab2 = calc_avcab(encab1,no), calc_avcab(encab2,no)
    # dH1 = H1 - <H1>, dH2 = H2-<H2>
    den1cab, den2cab = np.subtract(encab1,avencab1), np.subtract(encab2,avencab2)
    # dH1*dH2 
    d2encab = np.multiply(den1cab,den2cab)
    # <dH1*dH2>
    avd2encab = calc_avcab(d2encab,no)
    # <[dH1*dH2 - <dH1*dH2>]f(t)>
    #term1 = calc_avcab(np.multiply(d2encab,cab),nosplit*inputparam.segsplit)
    term1 = calc_avcab(np.multiply(d2encab,cab),no)
    #term2 = calc_avcab(np.multiply(avd2encab,cab),nosplit*inputparam.segsplit)
    term2 = calc_avcab(np.multiply(avd2encab,cab),no)
    av_d2cab = np.subtract(term1,term2)
    return av_d2cab

def calc_xderiv2(cab,encab1,encab2,no):
    """
    This calculates the second derivative as:
    f'' = <[dH^2 - <dH^2>]f>
    """
    #Second Derivative
    # <H1>, <H2>
    avencab1, avencab2 = calc_avcab(encab1,no), calc_avcab(encab2,no)
    # dH1 = H1 - <H1>, dH2 = H2-<H2>
    den1cab, den2cab = np.subtract(encab1,avencab1), np.subtract(encab2,avencab2)
    # dH1*dH2
    d2encab = np.multiply(den1cab,den2cab)
    term1 = calc_avcab(np.multiply(d2encab,cab),no)
    return term1

def second_order(cab,encab,no,item2):
    """
    This does the second derivative for all different options of item
    and blocking and error
    """
    print(item2)
    encab2 = np.multiply(energy[item2][:,None],no)
    d2 = calc_deriv2(cab,encab,encab2,no)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d2 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d2.append(calc_deriv2(cab[bstart:bend],encab[bstart:bend],encab2[bstart:bend],no[bstart:bend]))
        np.savetxt(subdir+"bl_"+str(block)+"_"+item+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d2[block]])
    err_d2 = Error(bl_d2)
    np.savetxt(subdir+item+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d2,err_d2])
    return

def second_cross(cab,encab,no):
    """
    This does the second derivative for all different options of item
    and blocking and error
    """
    for item2 in inp_n:
        print(item2)
        encab2 = np.multiply(energy[item2][:,None],no)
        d2 = calc_xderiv2(cab,encab,encab2,no)
        nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
        bl_d2 = []
        for block in range(inputparam.nblocks):
            bstart=block*nperb
            bend=(block+1)*nperb
            bl_d2.append(calc_xderiv2(cab[bstart:bend],encab[bstart:bend],encab2[bstart:bend],no[bstart:bend]))
            np.savetxt(subdir+"bl_"+str(block)+"_xd_"+item+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d2[block]])
        err_d2 = Error(bl_d2)
        np.savetxt(subdir+'xd_'+item+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d2,err_d2])
    return

def second_order_dist(cab,encab):
    """
    This does the second derivative for all different options of item
    and blocking and error
    """
    notmp = np.ones(np.shape(cab))
    for item2 in inp_n:
        print(item2)
        encab2=energy[item2]
        d2 = calc_deriv2(cab,encab,encab2,notmp)
        nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
        bl_d2 = []
        for block in range(inputparam.nblocks):
            bstart=block*nperb
            bend=(block+1)*nperb
            bl_d2.append(calc_deriv2(cab[bstart:bend],encab[bstart:bend],encab2[bstart:bend],notmp[bstart:bend]))
        err_d2 = Error(bl_d2)
        print("Second Derivative w.r.t. %s is: %13.9f +/- %13.9f" % (item2,d2,err_d2))
    return

def third_order(cab,encab,no):
    # Currently not supported.
    d3 = np.zeros(len(time))
    err_d3 = np.zeros(len(time))
    for item2 in inp_n:
        np.savetxt(subdir+item+'_'+item2+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d3,err_d3])

def fourth_order(cab,encab,no):
    # Currently not supported.
    d4 = np.zeros(len(time))
    err_d4 = np.zeros(len(time))
    for item2 in inp_n:
        np.savetxt(subdir+item+'_'+item2+'_'+item2+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d4,err_d4])

def Error(bl_corr):
    """
    Calcualtes the uncertainty in the right way from the blocks
    """
    return np.std(bl_corr,axis=0)*t_val



# Reads the command line arguments for the segments
parser = argparse.ArgumentParser()
parser.add_argument('-val', default=0, type=int, help='This is the number of the segment')
parser.add_argument('-opt', default=1, type=int, help='[1] init_segs [2] combine_segs [3] prew combine_segs')
parser.add_argument('-corr',default="c2", type=str,help='Correlation function of interest')
parser.add_argument('-mol', default="water",type=str,help="Molecule name")
parser.add_argument('-tnrm',default=-1, type=int, help="[-1] regular averaging, [1] special averaging")
parser.add_argument('-fname',default="flucts.inp", type=str, help="File that calls the types of energies")
parser.add_argument('-bin', default=0, type=int, help="[0] Non Binary [1] Binary")
parser.add_argument('-bfcorr', default="norm", type=str, help="norm_ is the default")
parser.add_argument('-time', default="real_time.dat", type=str, help="what file to read the time from")
parser.add_argument('-rehist', default=0, type=int, help="How to rebin the data")
parser.add_argument('-hmin', default=0, type=float, help="How to rebin the data")
parser.add_argument('-hmax', default=0, type=float, help="How to rebin the data")
parser.add_argument('-higher',default=1, type=int, help="[0] don't include higher derivatives")
parser.add_argument('-random', default=0, type=int, help="[0] don't randomize, [1] do randomize")
parser.add_argument('-weights', default="None", type=str, help="Name of correlation function for weights")
parser.add_argument('-wrule', default="linear", type=str, help="Options for weights, either 'linear' or 'square'")
parser.add_argument('-extra', default="None", type=str, help="Extra correlation function to read in")
args = parser.parse_args()
splitno     = args.val
option      = args.opt
corr_func   = args.corr
mol_name    = args.mol
tnrm        = args.tnrm
fname       = args.fname
binaryread  = args.bin
bfcorr      = args.bfcorr
tfile       = args.time
histbin     = args.rehist
hmin        = args.hmin
hmax        = args.hmax
randomize   = args.random
weights     = args.weights
higher      = args.higher
wrule       = args.wrule
extra_corr  = args.extra

# Read the input file
inputparam = user_input("input_file")

# Pull t-value from Student's T-Table
t_val = stats.t.ppf(0.975,inputparam.nblocks-1)/np.sqrt(inputparam.nblocks)

# Eventually these will be set by read_input
nosplit = int((inputparam.end_config-inputparam.start_config+1000)/(inputparam.segsplit*inputparam.sep_config))

# This calculates the beginning and end of the segment
fstart = splitno*inputparam.sep_config*inputparam.segsplit + inputparam.start_config
fend = (splitno+1)*inputparam.sep_config*inputparam.segsplit + inputparam.start_config

inp_n = np.genfromtxt(fname, usecols=0,dtype=str,unpack=True)
subdir="OUT/"

msd_flag = 0
if corr_func == "msdp":
    msd_flag = 1
    corr_func = "msd_py"



cab,no,wcab=[],[],[]
orig_cab,dist_av_cab = [],[]
encab=[]
extra_cab = []
if option == 1 and tnrm == -1:
    for i in range(fstart,fend,inputparam.sep_config):
        filename="FILES/"+str(i)+"/"+corr_func+"_"+str(i)+"_"+mol_name
        # This sets what values above tnrm to contribute, if -1, then disabled
        tmpcab = []
        # Note - support for binary read, much faster!
        #mj=open('missedjobs'+str(splitno), 'a')
        missedjobs=0
        if binaryread == 0: tmpcab = np.genfromtxt(filename+".dat", usecols=(1),dtype=(float),unpack=True)
        else: 
            try: 
                tmpcab = pickle.load(open(filename+".pckl",'rb'))
            except:
                missedjobs=1
                print(filename," has failed")
                #fdir = "FILES/"+str(i)
                #mj.write("mkdir "+fdir+"; cp in.nve "+fdir+"; cp nve.sh "+fdir+"; cd "+fdir+"; sed -i -e 's@AAA@"+str(i)+"@g' nve.sh; sbatch nve.sh; cd ../../\n")
        if missedjobs==1: exit()
        #mj.close()
        if histbin != 0:
            old = np.zeros(histbin)
            old[:np.array(tmpcab).shape[0]] = tmpcab
            tmpcab = old
        # Appends the corr func    
        cab.append(tmpcab)
    if not os.path.exists("TEMP"):
        os.makedirs("TEMP")
    if not os.path.exists("OUT"):
        os.makedirs("OUT")
    with open('TEMP/cab'+'_'+str(splitno)+'_'+corr_func+'_'+mol_name+'.pckl','wb') as g:
        pickle.dump(cab,g)
    # Only if using tnrm
    if tnrm != -1:
        with open('TEMP/no'+'_'+str(splitno)+'_'+corr_func+'_'+mol_name+'.pckl','wb') as h:
            pickle.dump(no,h)

elif option == 1 and tnrm != -1:
    for i in range(fstart,fend,inputparam.sep_config):
        filename="FILES/"+str(i)+"/"+corr_func+"_"+str(i)+"_"+mol_name
        # This sets what values above tnrm to contribute, if -1, then disabled
        tmpcab,tmpno = [],[]
        # Note - support for binary read - much faster!
        if binaryread == 0: tmpcab,tmpno = np.genfromtxt(filename+'.dat', usecols=(1,2),dtype=(float,float),unpack=True)
        else: 
            fname2 = "FILES/"+str(i)+"/"+bfcorr+"_"+str(i)+"_"+mol_name
            tmpcab = pickle.load(open(filename+".pckl",'rb'))
            tmpno  = pickle.load(open(fname2+".pckl", 'rb'))
        if histbin != 0:
            old = np.zeros(histbin)
            old[:np.array(tmpcab).shape[0]] = tmpcab
            tmpcab = old
        no.append((tmpno>tnrm)*1)
        # Appends the corr func
        cab.append(tmpcab)
    if not os.path.exists("TEMP"):
        os.makedirs("TEMP")
    if not os.path.exists("OUT"):
        os.makedirs("OUT")
    with open('TEMP/cab'+'_'+str(splitno)+'_'+corr_func+'_'+mol_name+'.pckl','wb') as g:
        pickle.dump(cab,g)
    # Only if using tnrm
    if tnrm != -1:
        with open('TEMP/no'+'_'+str(splitno)+'_'+corr_func+'_'+mol_name+'.pckl','wb') as h:
            pickle.dump(no,h)

elif option == 2:
    # Read in energies and put them in a dictionary
    energy={}
    print("Reading Energies")
    for item in inp_n:
        energy[item]=np.genfromtxt(item+'_init.out',usecols=0)
    # Read in the time, generated by real_time.py
    print("Reading Time")
    time = np.genfromtxt(tfile,usecols=0,unpack=True)
    # This loops over the segments and reads them in. 
    print("Unpickling Files")
    for i in range(nosplit):
        print(i,flush=True)
        tmpextra=[]
        tmpcab = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+corr_func+'_'+mol_name+'.pckl','rb'))
        if extra_corr != "None": tmpextra = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+extra_corr+'_'+mol_name+'.pckl','rb'))
        tmpweight=[]
        if weights != "None": tmpweight = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+weights+'_'+mol_name+'.pckl','rb'))
        if tnrm != -1: tmpno = pickle.load(open('TEMP/no'+'_'+str(i)+'_'+corr_func+'_'+mol_name+'.pckl','rb'))
        if i == 0:
            cab = tmpcab
            if extra_corr != "None": extra_cab = tmpextra
            if weights != "None": wcab = tmpweight
            if tnrm != -1: no  = tmpno
        else:
            cab=np.concatenate((cab,tmpcab),axis=0)
            if extra_corr != "None": extra_cab = np.concatenate((extra_cab,tmpextra),axis=0)
            if weights != "None": wcab = np.concatenate((wcab, tmpweight),axis=0)
            if tnrm != -1: no=np.concatenate((no,tmpno),axis=0)
    if msd_flag == 1:
        def linear(x,m,b):
            return m*x+b
        from scipy.optimize import curve_fit
        cab2 = []
        for msd in cab:
            cut = int(len(msd)*0.9)
            popt,pcov = curve_fit(linear, time[cut:], msd[cut:])
            m,b = popt
            cab2.append([m/.6])
        cab = cab2
        if len(np.shape(cab))!= 2: 
            np.savetxt("D_vals2.dat",np.c_[energy["e"],cab])
            cab=np.reshape(cab,(len(cab),1))
            if wrule=="inv": 
                cab = np.divide(1,cab)
                corr_func = corr_func + "inv"
            if wrule=="invlog":
                cab = np.log(np.divide(1,cab))
                corr_func = corr_func + "invlog"
            time=np.array([1])
    if "D" in corr_func:
        cab = np.reshape(cab, (len(cab),1))
        time = np.array([1])
    # Histograms if histbin is active
    # Slower!

    if histbin != 0:
        orig_cab = cab
        if wrule == "square": wcab = np.power(wcab,2)
        # Resets time if histogramming is active
        print("Binning Histogram",flush=True)
        split=(hmax-hmin)/float(histbin)
        time = []
        for b in range(histbin):
            time.append(b*split+split/2+hmin)
        # Resets cab if histogramming is active
        cab2 = []
        count = 0
        for tmpcab in cab:
            nonzero = tmpcab[tmpcab != 0]
            tmp, bedge= [],[]
            if count%1000 == 0: print("Reached histogram %d" % count, flush=True)
            if weights == "None": tmp,bedge = np.histogram(nonzero, bins=histbin,range=(hmin,hmax),density=True)
            else: tmp,bedge = np.histogram(nonzero, bins=histbin,range=(hmin,hmax),density=True, weights=wcab[count])
            # This next one is just a test
            #tmp = tmp / len(nonzero)
            cab2.append(tmp)
            count += 1
        cab = np.array(cab2,dtype=float)
    else:
        if extra_corr != "None":
            avextra = np.subtract(np.average(extra_cab,axis=1),np.average(extra_cab))
            cab = np.add(cab,np.power(avextra[:,np.newaxis],2))
        if wrule == "norm_zero":
            tmpav = np.average(cab,axis=0)
            cab = np.divide(cab,tmpav[0])
        elif wrule == "norm_reg":
            cab = np.divide(cab,cab[:,0,np.newaxis])
    # Sets Array for Normalization
    if tnrm == -1: no = np.ones(np.shape(cab))
    if histbin != 0: print("Histogramming Complete", flush=True)
    # Randomizes if option selected
    if randomize == 1:
        print("Randomizing the order")
        indices=np.arange(cab.shape[0])
        np.random.shuffle(indices)
        cab = cab[indices]
        no  = no[indices]
        for key in energy:
            energy[key]=energy[key][indices]

    if corr_func == "uh_theta": corr_func = "theta"
    if weights != "None": corr_func = corr_func + "_" + weights

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!THIS IS THE MAIN CORRELATION CALCULATION!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zero_order(cab,no)
    if histbin != 0:
        dist_av_cab = np.average(orig_cab,axis=1)
        zero_order_dist(dist_av_cab) 
    # Loops over energy items, reads them in to calculate correlations!
    for item in inp_n:
        print("Derivative for  %s" % item,flush=True)
        if histbin != 0: first_order_dist(dist_av_cab,energy[item])
        encab = np.multiply(energy[item][:,None],no)
        print(np.shape(encab))
        avcab = calc_avcab(cab,no)
        first_order(cab,encab,no)
        if higher == 1:
            second_order(cab,encab,no,item)
            second_cross(cab,encab,no)
            if histbin != 0: second_order_dist(dist_av_cab,energy[item])
            third_order(cab,encab,no)
            fourth_order(cab,encab,no)

elif option == 3:
    # For preweighted things already - first deriv only
    # Read in the time, generated by real_time.py
    print("Reading Time")
    time = np.genfromtxt(tfile,usecols=0,unpack=True)
    # This loops over the segments and reads them in.
    print("Unpickling Files")
    for i in range(nosplit):
        print(i,flush=True)
        tmpextra=[]
        tmpcab = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+corr_func+'_'+mol_name+'.pckl','rb'))
        if extra_corr != "None": tmpextra = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+extra_corr+'_'+mol_name+'.pckl','rb'))
        tmpweight=[]
        if weights != "None": tmpweight = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+weights+'_'+mol_name+'.pckl','rb'))
        if tnrm != -1: tmpno = pickle.load(open('TEMP/no'+'_'+str(i)+'_'+corr_func+'_'+mol_name+'.pckl','rb'))
        if i == 0:
            cab = tmpcab
            if extra_corr != "None": extra_cab = tmpextra
            if weights != "None": wcab = tmpweight
            if tnrm != -1: no  = tmpno
        else:
            cab=np.concatenate((cab,tmpcab),axis=0)
            if extra_corr != "None": extra_cab = np.concatenate((extra_cab,tmpextra),axis=0)
            if weights != "None": wcab = np.concatenate((wcab, tmpweight),axis=0)
            if tnrm != -1: no=np.concatenate((no,tmpno),axis=0)

    aven = 0.0 
    if weights != "None": aven = np.average(wcab[wcab!=0])
    if histbin != 0:
        orig_cab = cab
        if wrule == "square": wcab = np.power(wcab,2)
        # Resets time if histogramming is active
        print("Binning Histogram",flush=True)
        split=(hmax-hmin)/float(histbin)
        time = []
        for b in range(histbin):
            time.append(b*split+split/2+hmin)
        # Resets cab if histogramming is active
        cab2 = []
        origcab = []
        count = 0
        for tmpcab in cab:
            nonzero = tmpcab[tmpcab != 0]
            tmp, bedge= [],[]
            if count%1000 == 0: print("Reached histogram %d" % count, flush=True)
            orig,bedge = np.histogram(nonzero, bins=histbin,range=(hmin,hmax),density=False)
            nzwcab = wcab[count][tmpcab != 0]
            tmp,bedge = np.histogram(nonzero, bins=histbin,range=(hmin,hmax),density=False, weights=nzwcab)
            # This next one is just a test
            #tmp = tmp / len(nonzero)
            cab2.append(tmp)
            origcab.append(orig)
            count += 1
        cab = np.array(cab2,dtype=float)
        origcab = np.array(origcab,dtype=float)

        bl_cab,bl_dcab=[],[]
        nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
        for block in range(inputparam.nblocks):
            bstart=block*nperb
            bend=(block+1)*nperb
            blen = np.average(wcab[bstart:bend][wcab[bstart:bend]!=0],axis=0)
            orig = np.average(origcab[bstart:bend],axis=0)
            tmp = -np.subtract(np.average(cab[bstart:bend],axis=0),np.multiply(blen,orig))
            np.savetxt(subdir+"/bl_"+str(block)+"_"+mol_name+"_"+corr_func+".dat", np.c_[time,orig])
            np.savetxt(subdir+"/bl_"+str(block)+"_"+weights+"_"+mol_name+"_"+corr_func+".dat", np.c_[time,tmp])
            bl_dcab.append(tmp)
            bl_cab.append(orig)
        errcab = Error(bl_cab)
        errdcab = Error(bl_dcab)
        np.savetxt(subdir+"/"+mol_name+"_"+corr_func+".dat", np.c_[time,np.average(origcab,axis=0),errcab])
        fincab = -np.subtract(np.average(cab,axis=0),np.multiply(aven,np.average(origcab,axis=0)))
        np.savetxt(subdir+"/"+weights+"_"+mol_name+"_"+corr_func+".dat", np.c_[time,fincab,errdcab])
    else:
        if extra_corr != "None":
            avextra = np.subtract(np.average(extra_cab,axis=1),np.average(extra_cab))
            cab = np.add(cab,np.power(avextra[:,np.newaxis],2))
        if wrule == "norm_zero":
            tmpav = np.average(cab,axis=0)
            cab = np.divide(cab,tmpav[0])
        elif wrule == "norm_reg":
            cab = np.divide(cab,cab[:,0,np.newaxis])
    
    
