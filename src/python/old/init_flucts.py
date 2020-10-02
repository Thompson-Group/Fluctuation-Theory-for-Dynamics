import numpy as np
import sys
from read_input import input

def BLOCK_ENERGY(energy,start,end,itemindex):
    eav = 0
    for i in range(start*sep, end*sep):
        eav += energy[itemindex][i]
    eav = eav/(end*sep-start*sep)
    return eav

inputparam = input("input_file")
sep = 500
num_segs = int(inputparam.num_files/float(sep))
segs_per_block = np.floor(num_segs/float(inputparam.nblocks))

val = sys.argv[1]
fname = str(sys.argv[2])
corr_name = str(sys.argv[3])
mol_name = str(sys.argv[4])

# Read In Time
time = np.genfromtxt('real_time.dat', usecols=0)
time = [0] + time
time = time[:inputparam.num_times]

# Read In Fluctuations
inp_n, inp_c=np.genfromtxt(fname, usecols=(0,1), dtype=(str,int),unpack=True)

# Initialize Arrays
d1=[sep]
d2=[sep]
d3=[sep]
d4=[sep]
corr=[]
d1corr=[]
d2corr=[]
d3corr=[]
d4corr=[]

for i in range(inputparam.num_times):
    corr.append(0)
    d1corr.append(0)
    d2corr.append(0)
    d3corr.append(0)
    d4corr.append(0)

energy=[]
tcab=[]

print("The Argument Provided is %s:" % val)
if val == '-h':
    print("This function needs a single integer input to run.")
    print("     Usage: python init_average.py sep fname corr_func")
    exit()
else:
    itmp=int(val)
    print itmp
    inputparam = input("input_file")
    fstart = itmp*sep
    fend = (itmp+1)*sep
    print("Starting Trajectory: %s\nEnding Trajectory: %s" % (fstart,fend))
    # Read in Correlation Functions
    fnames = np.genfromtxt('file_names', dtype='string', unpack=True)
    print("Choosing %s of %s files" % (sep, len(fnames)))
    fcab=["FILES/"+str(s)+"/"+str(corr_name)+"_"+str(s)+"_"+str(mol_name)+".dat" for s in fnames[fstart:fend]]
    item1count=0
    for item1 in inp_n:
        energy.append(np.genfromtxt(item1+'_init.out'))
    for i in range(sep):
        tcab.append(np.genfromtxt(fcab[i],usecols=1,unpack=True))
    print len(tcab)
    for item1 in inp_n:
        print item1
        item2count=0
        for item2 in inp_n:
            corr=np.zeros((inputparam.num_times))
            bl_corr = np.zeros((inputparam.num_times))
            bl_w1corr  = np.zeros((inputparam.num_times))
            bl_w2corr  = np.zeros((inputparam.num_times))
            bl_w3corr  = np.zeros((inputparam.num_times))
            bl_w4corr  = np.zeros((inputparam.num_times))
            w1corr=np.zeros((inputparam.num_times))
            w2corr=np.zeros((inputparam.num_times))
            w3corr=np.zeros((inputparam.num_times))
            w4corr=np.zeros((inputparam.num_times))
            for block in range(inputparam.nblocks):
                bstart = int(block*segs_per_block)
                bend = int((block+1)*segs_per_block)
                bdist = bend - bstart
                bl_e1av = BLOCK_ENERGY(energy,bstart,bend,item1count)
                bl_e2av = BLOCK_ENERGY(energy,bstart,bend,item2count)
                for i in range(sep):
                    enum = fstart + i
                    d1=energy[item1count][enum]-bl_e1av
                    d2=energy[item2count][enum]-bl_e2av
                    bl_corr = bl_corr + tcab[i]
                    bl_w1corr=bl_w1corr+tcab[i]*d1
                    bl_w2corr=bl_w2corr+tcab[i]*d1*d2
                    bl_w3corr=bl_w3corr+tcab[i]*d1*d2*d2
                    bl_w4corr=bl_w4corr+tcab[i]*d1*d2*d2*d2 
                bl_corr[:] = [x / float(sep) for x in bl_corr]
                bl_w1corr[:] = [x / float(sep) for x in bl_w1corr]
                bl_w2corr[:] = [x / float(sep) for x in bl_w2corr]
                bl_w3corr[:] = [x / float(sep) for x in bl_w3corr]
                bl_w4corr[:] = [x / float(sep) for x in bl_w4corr] 
                np.savetxt('SEG/bl_'+str(block)+'_seg_'+str(int(val))+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', np.c_[time, bl_corr, bl_w1corr, bl_w2corr, bl_w3corr, bl_w4corr])

            # Calculate Average Energy
            seg_start = 0
            seg_end = num_segs
            seg_dist = seg_end - seg_start
            e1av = BLOCK_ENERGY(energy,seg_start,seg_end,item1count)
            e2av = BLOCK_ENERGY(energy,seg_start,seg_end,item2count)
            for i in range(sep):
                enum=fstart+i
                corr=corr+tcab[i]
                d1=energy[item1count][enum]-e1av
                d2=energy[item2count][enum]-e2av
                w1corr=w1corr+tcab[i]*d1
                w2corr=w2corr+tcab[i]*d1*d2
                w3corr=w3corr+tcab[i]*d1*d2*d2
                w4corr=w4corr+tcab[i]*d1*d2*d2*d2
            corr[:] = [x / float(sep) for x in corr]
            w1corr[:] = [x / float(sep) for x in w1corr]
            w2corr[:] = [x / float(sep) for x in w2corr]
            w3corr[:] = [x / float(sep) for x in w3corr]
            w4corr[:] = [x / float(sep) for x in w4corr]
            np.savetxt('SEG/seg_'+str(int(val))+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', np.c_[time, corr, w1corr, w2corr, w3corr, w4corr])
            item2count+=1
        item1count+=1

