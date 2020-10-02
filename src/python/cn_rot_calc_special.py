#!/usr/bin/env python
import numpy as np
import sys
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-atom1', default=1,type=int,help="At what timestep should jumps be split into early or late.")
parser.add_argument('-atom2', default=2,type=int,help="At what timestep should jumps be split into early or late.")
args = parser.parse_args()

atom1, atom2 = args.atom1, args.atom2

atomselect=np.array([atom1,atom2])


def read_inp(inpfile):
    """
    This reads the input file that is used by the direct method
    for the calculations that define the initial parameters.
    """
    inp = open(inpfile, 'r')
    params = {}
    lines = []
    for line in inp:
        lines.append(line)
    params["nfile"]  = str(lines[1]).strip()
    params["ntimes"] = int(lines[3].split()[0])
    params["dt"]     =float(lines[3].split()[1])
    params["volume"] = float(lines[5])
    params["mol_name"] = str(lines[7]).strip()
    inp.close()

    params["L"] = params["volume"]**(1/3.)
    print("Input parameters read in")
    for key in params:
        print("%8s = %8s" % (key, str(params[key])))
    return params

def dist_components(L,rA,rB):
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

def read_frame(f,atomselect,params):
    """
    Reads a single frame at a time for an inputfile (already opened as f)
    """
    r_atom = []
    count,index,nmols = 0,0,0
    while True:
        line = f.readline()
        if len(line.split()) == 1 and count > 0:
            r_atom = np.array(r_atom)
            tmpe,tmpr=dist_components(params["L"],r_atom,r_atom)
            return tmpr, -tmpe
        if len(line.split()) == 4:
            count += 1
            if count in atomselect:
                r_atom.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
    sys.exit("Error: File read was unsucessful. - never read the right type of file")

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

def time_loop(params,filename,atomselect):
    """
    The easiest cn code ever, since it doesn't need an average over molecules.
    """
    f = open(filename,'r')
    r, e, eo = [],[],[]
    c1,c2,c3 = np.zeros(params["ntimes"]), np.zeros(params["ntimes"]), np.zeros(params["ntimes"])
    for t in range(params["ntimes"]):
        r, e = read_frame(f, atomselect,params)
        if t == 0:
            eo = e[0][1]
        e=e[0][1]
        c1[t], c2[t], c3[t] = calc_cn(eo,e)
    f.close()
    return c1,c2,c3



if __name__ == "__main__":
    params = read_inp("corr_calc.in")
    filename = 'traj_'+str(params["nfile"])+'_'+str(params["mol_name"])+'.xyz'
    c1,c2,c3 = time_loop(params, filename, atomselect)
    steps = np.linspace(0,params["ntimes"]-1,num=params["ntimes"])
    np.savetxt('spec_c1_'+str(params["nfile"])+'_'+str(params["mol_name"])+'.dat', np.c_[steps,c1], fmt="%2.5f")
    np.savetxt('spec_c2_'+str(params["nfile"])+'_'+str(params["mol_name"])+'.dat', np.c_[steps,c2], fmt="%2.5f")
    np.savetxt('spec_c3_'+str(params["nfile"])+'_'+str(params["mol_name"])+'.dat', np.c_[steps,c3], fmt="%2.5f")
    

