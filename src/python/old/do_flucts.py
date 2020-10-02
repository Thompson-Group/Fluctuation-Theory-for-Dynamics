import numpy as np
import sys
from scipy import stats
from read_input import input

"""
This is a python program to take the weighted correlation functions calculated in init_flucts.py and turn them into the derivative correlation functions.
"""

def BLOCK_ENERGY(energy,start,end,itemindex):
    eav = 0
    for i in range(start*sep, end*sep):
        eav += energy[itemindex][i]
    eav = eav/(end*sep-start*sep)
    return eav

def BLOCK_DENERGY(energy,start,end,itemindex,eav,n):
    deav = 0
    for i in range(start*sep, end*sep):
        deav += (energy[itemindex][i]-eav)**n
    deav = deav/(end*sep-start*sep)
    return deav

def BLOCK_ARRAY_ERR(barray, noblocks):
    sigma = []
    for i in range(len(barray)):
        sigma.append(np.std(barray[i])*t_val)
    return sigma

def NORM(array, n):
    return array/float(n)

"""
Derivative for the first derivative
C'(t) = -<dh(0)A(0)B(t)>
"""
def FIRST_DERIV(w1corr):
    d1corr = -w1corr
    return d1corr

"""
Derivative for the second derivative
C''(t) = <dh(0)^2A(0)B(t)>-<dh^2>C(t)
"""
def SECOND_DERIV(corr, w2corr, d2av):
    d2corr = w2corr - d2av*corr
    return d2corr

"""
Derivative for the third derivative
C'''(t) = -<dh(0)^3A(0)B(t)> + <dh^3>C(t) - 3<dh^3>C'(t)
"""
def THIRD_DERIV(corr, d1corr, w3corr, d2av, d3av):
    d3corr = -w3corr - 3*d2av*d1corr + d3av*corr
    return d3corr
"""
Derivative for the fourth derivative
C''''(t) = <[dh(0)^4 -<dh^4>]A(0)B(t)> - 6*<dh^2>C''(t)+4<dh^3>C'(t)
"""
def FOURTH_DERIV(corr, d1corr, d2corr, w4corr, d2av, d3av, d4av):
    d4corr = w4corr - 6*d2av*d2corr + 4*d3av*d1corr - d4av*corr
    return d4corr

def RATIO(corr, d1corr):
    ea = np.divide(d1corr[1:],corr[1:])
    ea=np.insert(ea,0,0.0)
    return ea
    


# Import Input File Parameters
inputparam = input('input_file')

# Read in command line args
if len(sys.argv) != 4:
    print("Usage: python do_flucts.py fname corr_name mol_name")
    exit(1)
fname = str(sys.argv[1])
corr_name = str(sys.argv[2])
mol_name = str(sys.argv[3])

# Pull t-value from Student's T-Table
t_val = stats.t.ppf(0.975,inputparam.nblocks-1)/np.sqrt(inputparam.nblocks)

#Calculate number of segments
sep = 500
num_segs = int(inputparam.num_files/float(sep))
segs_per_block = np.floor(num_segs/float(inputparam.nblocks))
print("There are %s total segments" % num_segs)
print("There are %s blocks" % inputparam.nblocks)
print("There are %s segs_per_block" % segs_per_block)

# Read in corr-names
inp_n, inp_c = np.genfromtxt(fname, usecols=(0,1), dtype=(str,int), unpack=True)

# Initialize Arrays
time   = np.zeros(inputparam.num_times)
blcorr   = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, inputparam.num_times))
blw1corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, inputparam.num_times))
blw2corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, inputparam.num_times))
blw3corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, inputparam.num_times))
blw4corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, inputparam.num_times))
corr   = np.zeros((len(inp_n),len(inp_n), num_segs, inputparam.num_times))
w1corr = np.zeros((len(inp_n),len(inp_n), num_segs, inputparam.num_times))
w2corr = np.zeros((len(inp_n),len(inp_n), num_segs, inputparam.num_times))
w3corr = np.zeros((len(inp_n),len(inp_n), num_segs, inputparam.num_times))
w4corr = np.zeros((len(inp_n),len(inp_n), num_segs, inputparam.num_times))

energy = []

# Read Segments In
item1count=0
print("---Start Read---")
for item1 in inp_n:
    print("The current Item is: %s" % item1)
    item2count=0
    for item2 in inp_n:
        for seg in range(num_segs):
            # 4D Arrays of shape [item1][item2][seg][time]
            time,corr[item1count][item2count][seg],w1corr[item1count][item2count][seg],w2corr[item1count][item2count][seg],w3corr[item1count][item2count][seg],w4corr[item1count][item2count][seg]=np.genfromtxt('SEG/seg_'+str(seg)+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', usecols=(0,1,2,3,4,5), unpack=True)
            for block in range(inputparam.nblocks):
                blcorr[block][item1count][item2count][seg],blw1corr[block][item1count][item2count][seg],blw2corr[block][item1count][item2count][seg],blw3corr[block][item1count][item2count][seg],blw4corr[block][item1count][item2count][seg] = np.genfromtxt('SEG/bl_'+str(block)+'_seg_'+str(seg)+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', usecols=(1,2,3,4,5), unpack=True)

        item2count += 1
    item1count += 1
print("---End Read---")
print("---Start Blocking---")
for item1 in inp_n:
    energy.append(np.genfromtxt(item1+'_init.out'))
item1count = 0
for item1 in inp_n:
    print("Starting loop for %s" % item1)
    item2count = 0
    for item2 in inp_n:
        print("Calculating %s %s" % (item1,item2))
        item3 = item2
        item4 = item2
        # Zero Block Weighted Arrays
        bl_w1corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_w2corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_w3corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_w4corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        # Zero Block Derivative Arrays 
        bl_corr    = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_d1corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_d2corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_d3corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_d4corr  = np.zeros((inputparam.nblocks, inputparam.num_times))
        bl_ea      = np.zeros((inputparam.nblocks, inputparam.num_times))
        # Zero Tot Weighted Arrays
        tot_w1corr = np.zeros(inputparam.num_times)
        tot_w2corr = np.zeros(inputparam.num_times)
        tot_w3corr = np.zeros(inputparam.num_times)
        tot_w4corr = np.zeros(inputparam.num_times)
        # Zero Tot Derivative Arrays
        tot_corr   = np.zeros(inputparam.num_times)
        tot_d1corr = np.zeros(inputparam.num_times)
        tot_d2corr = np.zeros(inputparam.num_times)
        tot_d3corr = np.zeros(inputparam.num_times)
        tot_d4corr = np.zeros(inputparam.num_times)
        tot_ea     = np.zeros(inputparam.num_times)
        # Zero Err Weighted Arrays
        err_w1corr = np.zeros(inputparam.num_times)
        err_w2corr = np.zeros(inputparam.num_times)
        err_w3corr = np.zeros(inputparam.num_times)
        err_w4corr = np.zeros(inputparam.num_times)
        # Zero Err Derivative Arrays
        err_corr   = np.zeros(inputparam.num_times)
        err_d1corr = np.zeros(inputparam.num_times)
        err_d2corr = np.zeros(inputparam.num_times)
        err_d3corr = np.zeros(inputparam.num_times)
        err_d4corr = np.zeros(inputparam.num_times)
        err_ea     = np.zeros(inputparam.num_times)
        # Sum Segments into blocks
        for block in range(inputparam.nblocks):
            print("     BLOCK %s" % block)
            # Block indices calculation
            bstart = int(block*segs_per_block)
            bend   = int((block+1)*segs_per_block)
            bdist  = bend - bstart
            # Calculate Average Flucts
            e1_av  =   BLOCK_ENERGY(energy, bstart, bend, item1count)
            d1_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,1)
            d2_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,2)
            d3_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,3)
            d4_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,4)

            for i in range(inputparam.num_times):
                for seg in range(bstart,bend):
                    # Need to normalize.
                    bl_corr[block][i]   +=   blcorr[block][item1count][item2count][seg][i]
                    bl_w1corr[block][i] += blw1corr[block][item1count][item2count][seg][i]
                    bl_w2corr[block][i] += blw2corr[block][item1count][item2count][seg][i] 
                    bl_w3corr[block][i] += blw3corr[block][item1count][item2count][seg][i]
                    bl_w4corr[block][i] += blw4corr[block][item1count][item2count][seg][i]
                # Average over segments
                bl_corr[block][i]   = NORM(bl_corr[block][i], bdist)
                bl_w1corr[block][i] = NORM(bl_w1corr[block][i], bdist)
                bl_w2corr[block][i] = NORM(bl_w2corr[block][i], bdist)
                bl_w3corr[block][i] = NORM(bl_w3corr[block][i], bdist)
                bl_w4corr[block][i] = NORM(bl_w4corr[block][i], bdist)
                # Calculate Derivatives
                bl_d1corr[block][i] = FIRST_DERIV(bl_w1corr[block][i])
                bl_d2corr[block][i] = SECOND_DERIV(bl_corr[block][i],bl_w2corr[block][i],d2_av)
                bl_d3corr[block][i] = THIRD_DERIV(bl_corr[block][i], bl_d1corr[block][i],bl_w3corr[block][i],d2_av, d3_av)
                bl_d4corr[block][i] = FOURTH_DERIV(bl_corr[block][i], bl_d1corr[block][i], bl_d2corr[block][i], bl_w4corr[block][i], d2_av, d3_av, d4_av)
            # Calculate ratio function
            bl_ea[block] = RATIO(bl_corr[block], bl_d1corr[block])
            # Sets File Names
            bl_name   = "bl_"+str(block)+"_"+mol_name+"_"+corr_name+".dat"
            bl_d1name = "bl_"+str(block)+"_"+item1+"_"+mol_name+"_"+corr_name+".dat"
            bl_d2name = "bl_"+str(block)+"_"+item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat"
            bl_d3name = "bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+mol_name+'_'+corr_name+".dat"
            bl_d4name = "bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+'_'+corr_name+".dat"
            # Writes to file
            np.savetxt(bl_name,   np.c_[time, bl_corr[block], bl_ea[block]],      fmt='%s')
            np.savetxt(bl_d1name, np.c_[time, bl_d1corr[block],bl_w1corr[block]], fmt='%s')
            np.savetxt(bl_d2name, np.c_[time, bl_d2corr[block],bl_w2corr[block]], fmt='%s')
            np.savetxt(bl_d3name, np.c_[time, bl_d3corr[block],bl_w3corr[block]], fmt='%s')
            np.savetxt(bl_d4name, np.c_[time, bl_d4corr[block],bl_w4corr[block]], fmt='%s')
        # Calculate Uncertainties
        err_corr    =   np.array(bl_corr).std(0)
        err_corr    =   [x * t_val for x in err_corr]
        err_w1corr  =   np.array(bl_w1corr).std(0)
        err_w1corr  =   [x * t_val for x in err_w1corr]
        err_w2corr  =   np.array(bl_w2corr).std(0)
        err_w2corr  =   [x * t_val for x in err_w2corr]
        err_w3corr  =   np.array(bl_w3corr).std(0)
        err_w3corr  =   [x * t_val for x in err_w3corr]
        err_w4corr  =   np.array(bl_w4corr).std(0)
        err_w4corr  =   [x * t_val for x in err_w4corr]
        err_d1corr  =   np.array(bl_d1corr).std(0)
        err_d1corr  =   [x * t_val for x in err_d1corr]
        err_d2corr  =   np.array(bl_d2corr).std(0)
        err_d2corr  =   [x * t_val for x in err_d2corr]
        err_d3corr  =   np.array(bl_d3corr).std(0)
        err_d3corr  =   [x * t_val for x in err_d3corr]
        err_d4corr  =   np.array(bl_d4corr).std(0)
        err_d4corr  =   [x * t_val for x in err_d4corr]
        err_ea      =   np.array(bl_ea).std(0)
        err_ea      =   [x * t_val for x in err_corr]
        seg_start = 0
        seg_end   = num_segs
        seg_dist  = seg_end - seg_start
        e1_av = BLOCK_ENERGY(energy, seg_start, seg_end, item1count)
        d1_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 1)
        d2_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 2)
        d3_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 3)
        d4_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 4)
        print("d1_av = %s" % d1_av)
        print("d2_av = %s" % d2_av)
        print("d3_av = %s" % d3_av)
        print("d4_av = %s" % d4_av)
        for i in range(inputparam.num_times):
            for seg in range(num_segs):
                tot_corr[i]   += corr[item1count][item2count][seg][i]
                tot_w1corr[i] += w1corr[item1count][item2count][seg][i]
                tot_w2corr[i] += w2corr[item1count][item2count][seg][i]
                tot_w3corr[i] += w3corr[item1count][item2count][seg][i]
                tot_w4corr[i] += w4corr[item1count][item2count][seg][i]
            # Normalize
            tot_corr[i] = NORM(tot_corr[i], seg_dist) 
            tot_w1corr[i] = NORM(tot_w1corr[i], seg_dist)
            tot_w2corr[i] = NORM(tot_w2corr[i], seg_dist)
            tot_w3corr[i] = NORM(tot_w3corr[i], seg_dist)
            tot_w4corr[i] = NORM(tot_w4corr[i], seg_dist)
            # Calculate Derivatives
            tot_d1corr[i] = FIRST_DERIV(tot_w1corr[i])
            tot_d2corr[i] = SECOND_DERIV(tot_corr[i],tot_w2corr[i],d2_av)
            tot_d3corr[i] = THIRD_DERIV(tot_corr[i], tot_d1corr[i],tot_w3corr[i],d2_av, d3_av)
            tot_d4corr[i] = FOURTH_DERIV(tot_corr[i], tot_d1corr[i], tot_d2corr[i], tot_w4corr[i], d2_av, d3_av, d4_av)
        tot_ea = RATIO(tot_corr,tot_d1corr)
        tot_name   = mol_name+"_"+corr_name+".dat"
        tot_d1name = item1+"_"+mol_name+"_"+corr_name+".dat"
        tot_d2name = item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat"
        tot_d3name = item1+"_"+item2+"_"+item3+"_"+mol_name+"_"+corr_name+".dat"
        tot_d4name = item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+"_"+corr_name+".dat"
        np.savetxt(tot_name,   np.c_[time, tot_corr, err_corr, tot_ea, err_ea], fmt='%s')
        np.savetxt(tot_d1name, np.c_[time, tot_d1corr, err_d1corr, tot_w1corr, err_w1corr], fmt='%s')
        np.savetxt(tot_d2name, np.c_[time, tot_d2corr, err_d2corr, tot_w2corr, err_w2corr], fmt='%s')
        np.savetxt(tot_d3name, np.c_[time, tot_d3corr, err_d3corr, tot_w3corr, err_w3corr], fmt='%s')
        np.savetxt(tot_d4name, np.c_[time, tot_d4corr, err_d4corr, tot_w4corr, err_w4corr], fmt='%s')
        
        item2count+=1
    item1count+=1


            
            




               







