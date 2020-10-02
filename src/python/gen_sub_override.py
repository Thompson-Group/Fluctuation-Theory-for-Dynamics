# This python script generates the needed sub.sh, nve.sh, and job_array.sh scripts in the simulation directory i.e. /sim/
"""
This python script generates the needed files in the simulation directory.
Requires an input_file
This produces the following output:
    run_array.sh
    nve.sh
    init_segments.sh
    combine_segments.sh
"""

# These set global variables and should be changed based on your system
lammpsload="module load compiler/intel/18 intel-mpi/18\nmodule load lammps/3Mar2020"
lammpsruncmd="mpirun lmp_mpi < in.nve -screen none"
cp2kload="module load cp2k/6.0/popt"
cp2kruncmd="mpirun cp2k.popt in.nve.cp2k"
partition="sixhour"
#cp2kruncmd="    srun -n 128 -c 8 --cpu_bind=cores  /global/homes/c/cjmundy/SOURCEFORGE_CP2K/cp2k-5.1/cp2k/exe/CRAY-cori-Haswell-gnu/cp2k.psmp in.nve.cp2k > cp2k.out"


def read_machine():
    f = open('machine.name','r')
    machine = str(f.readline().strip())
    path = str(f.readline().strip())
    return machine,path

def read_override():
    f=open(".flag_override",'r')
    fname=f.readline().strip()
    f.close()
    override = open(fname,'r')
    lines = override.readlines()
    return lines


import numpy as np
import sys
from src.python.read_input import user_input
import os


print("WARNING: Override Enabled, please make appropriate changes in the file specified")
override=read_override()

machine, homepath = read_machine()

lines=[]
with open(homepath+"/src/dependencies/"+machine+"_header.dat",'r') as mfile:
    for line in mfile:
        lines.append(line)

# Calls the read_input class
inputparam = user_input("input_file")

njobs = int(np.ceil(inputparam.num_files/float(inputparam.num_rpj)))
# Generates a single job array that runs 50 jobs each.
dcn_file = "run_array.sh"
dcn = open(dcn_file, 'w')
for item in lines:
    dcn.write(item)
dcn.write("#SBATCH --array 0-%s\n" % (njobs-1))
dcn.write("\n")
dcn.write("module load Dir_Calc_Fluct\n")
dcn.write("cd $SLURM_SUBMIT_DIR\n")
if inputparam.prog == "LAMMPS":
    dcn.write("%s\n" % lammpsload)
elif inputparam.prog == "CP2K":
    dcn.write("%s\n" % cp2kload)
else:
    dcn.write("Error: Incorrect program in input_file\n")
    exit(1)

if machine == "CORI":
    dcn.write("export OMP_NUM_THREADS=8\n")
    dcn.write("export OMP_PLACES=threads\n")
    dcn.write("export OMP_PROC_BIND=spread\n")

dcn.write("\n")
dcn.write("SEP=%s\n" % int(inputparam.sep_config))
dcn.write("START=%s\n" % int(inputparam.start_config))
dcn.write("\n")
dcn.write("START_JOBS=1\n")
dcn.write("NUM_RPJ=%s\n" % inputparam.num_rpj)
dcn.write("\n")
dcn.write("for (( N = $START_JOBS; N <= $NUM_RPJ; N++ ))\n")
dcn.write("do\n")
dcn.write("    CUR=$(( START + N*SEP + SEP*NUM_RPJ*SLURM_ARRAY_TASK_ID - SEP))\n")
dcn.write("    echo $CUR\n")
dcn.write("    mkdir -p FILES/$CUR\n")
if inputparam.prog == "LAMMPS":
    dcn.write("    cp in.nve nve.sh FILES/$CUR/\n")
elif inputparam.prog == "CP2K":
    dcn.write("    cp in.nve.cp2k nve.sh FILES/$CUR/\n")
else:
    dcn.write("Error: Incorrect program in input_file")
dcn.write("    cd FILES/$CUR\n")
dcn.write("    \n\n")
dcn.write("    sed -i -e 's@direct_calc_nve@nve_'$CUR'@g' nve.sh\n")
if inputparam.prog == "LAMMPS":
    dcn.write("    sed -i -e 's@restart.file@../../RESTART/restart.'$CUR'@g' in.nve\n")
    for j in range(0,inputparam.num_molecs):
        dcn.write("    sed -i -e 's@traj.file%s@traj_'$CUR'_%s.xyz@g' in.nve\n" % (str(j),inputparam.molec[j]))
    dcn.write("    %s\n" % lammpsruncmd)
elif inputparam.prog == "CP2K":
    if inputparam.cab == "IONPAIRING":
        dcn.write("    vel_reselect.py 6.94 18.998 298.15 restart.$CUR\n")
        dcn.write("    sed -i -e 's@vel.file@vel_'$CUR'_%s.vxyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    dcn.write("    sed -i -e 's@restart.file@RESTART/restart.'$CUR'@g' in.nve.cp2k\n")
    dcn.write("    sed -i -e 's@traj.file0@traj_'$CUR'_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    dcn.write("    %s\n" % runcmd)

dcn.write("    \n\n")
for line in override:
    dcn.write("     %s" % line)
dcn.write('    touch corr_complete\n')
dcn.write("    find . -type f -not \( -name '*pckl' -or -name '*lammps' \) -delete\n")
dcn.write("    cd ../../\n")
dcn.write("done\n")
dcn.close()


# Generate NVE Input File
#   This is the file that is generated in every directory that needs it if
#   something doesn't run correctly.
nve_file="nve.sh"
nve=open(nve_file, 'w')

for item in lines:
    if "output" in line:
        nve.write("#SBATCH --output=nve.out")
    else:
        nve.write(item)
nve.write("module load Dir_Calc_Fluct\n")
if inputparam.prog == "LAMMPS":
    nve.write('%s\n\n' % lammpsload)
elif inputparam.prog == "CP2K":
    nve.write('%s\n\n' % cp2kload)
if machine == "CORI":
    nve.write("export OMP_NUM_THREADS=8\n")
    nve.write("export OMP_PLACES=threads\n")
    nve.write("export OMP_PROC_BIND=spread\n")
nve.write('echo Time is `date` > array.o\n')
nve.write('echo Directory is `pwd` >> array.o\n\n\n')

if inputparam.prog == "LAMMPS":
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve\n")
    for j in range(0,inputparam.num_molecs):
        nve.write("sed -i -e 's@traj.file%s@traj_AAA_%s.xyz@g' in.nve\n" % (str(j),inputparam.molec[j]))
    nve.write('%s\n\n\n' % lammpsruncmd)
elif inputparam.prog == "CP2K":
    if inputparam.cab == "IONPAIRING":
        nve.write("vel_reselect.py 6.94 18.998 298.15 restart.$CUR\n")
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve.cp2k\n")
    nve.write("sed -i -e 's@traj.file0@traj_AAA_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    nve.write("%s\n" % cp2kruncmd)

for line in override:
    nve.write("     %s" % line)

nve.write('touch corr_complete\n')
nve.write('echo Ending Time is `date` >> array.o\n')


nve.close()

# Gen Init Array
init_file="init_segments.sh"
iarr = open(init_file,'w')
sep = int(inputparam.num_files/inputparam.segsplit)
if inputparam.num_files % inputparam.nblocks != 0:
    print("Num files and Nblocks do not divide evenly, aborting")
    sys.exit(1)
iarr.write("#!/bin/bash\n")
iarr.write("#SBATCH --job-name=init_array\n")
iarr.write("#SBATCH --partition=%s\n" % partition)
iarr.write("#SBATCH --output=init_array.out\n")
iarr.write("#SBATCH --nodes=1\n")
iarr.write("#SBATCH --ntasks-per-node=2\n")
iarr.write("#SBATCH --constraint=intel\n")
iarr.write("#SBATCH --mem=30G\n")
iarr.write("#SBATCH --time=06:00:00\n")
iarr.write("#SBATCH --array=0-%s\n" % (sep-1))

iarr.write("module load Dir_Calc_Fluct\n")

iarr.write("cd $SLURM_SUBMIT_DIR\n")
iarr.write(" parse_and_combine.py -opt 1 -bin 1 -val $SLURM_ARRAY_TASK_ID -fname flucts.inp -corr msd -mol %s \n" % inputparam.molec[0])

iarr.close()
print("WARNING: Because override is enabled, init_segments.sh will need to be modified to include desired commands")

# Gen Combine Segments Array
# This submits the combine segmentspython code which does the final analysis.

do_file="combine_segments.sh"
darr = open(do_file, 'w')
darr.write("#!/bin/bash\n")
darr.write("#SBATCH --job-name=do_flucts\n")
darr.write("#SBATCH --partition=%s\n" % partition)
darr.write("#SBATCH --output=do_flucts.out\n")
darr.write("#SBATCH --nodes=1\n")
darr.write("#SBATCH --ntasks-per-node=20\n")
darr.write("#SBATCH --constraint=intel\n")
darr.write("#SBATCH --mem=100G\n")
darr.write("#SBATCH --time=06:00:00\n")
darr.write("module load Dir_Calc_Fluct\n")
darr.write("srun -N1 -n1 -c1 --mem=20G --exclusive parse_and_combine.py -opt 2 -fname flucts.inp -corr msd -mol %s &\n" % (inputparam.molec[0]))
darr.write("wait\n")
darr.close()
print("WARNING: Because override is enabled, combine_segments.sh will need to be modified to include desired commands. An example script, with an example command has been pre-populated for simplicity.")

# Generate Flucts Script

import os

ival,icol=np.genfromtxt("flucts.inp", usecols=(0,1),dtype=(str,int), unpack=True)
count = 0
for item in ival:
    if "D" not in item and "A" not in item:
        count+=1


flc=open('grabfluctsub.sh','w')
flc.write("#!/bin/bash\n")
flc.write("#SBATCH --job-name=grabflucts_array\n")
flc.write("#SBATCH --partition=sixhour\n")
flc.write("#SBATCH --output=grabflucts_array%A_%a.out\n")
flc.write("#SBATCH --nodes=1\n")
flc.write("#SBATCH --ntasks-per-node=2\n")
flc.write("#SBATCH --constraint=intel\n")
flc.write("#SBATCH --mem=30G\n")
flc.write("#SBATCH --time=06:00:00\n")
flc.write("#SBATCH --array=0-%d\n" % int(count-1))

flc.write("module load Dir_Calc_Fluct\n")

flc.write("cd $SLURM_SUBMIT_DIR\n")

flc.write("echo Time is `date` > array_$SLURM_ARRAY_TASK_ID.o\n")
flc.write("echo Directory is `pwd` >> array_$SLURM_ARRAY_TASK_ID.o\n")

flc.write("grab_flucts.py flucts.inp $SLURM_ARRAY_TASK_ID\n")

flc.close()

pwt=open('comboflucts.sh', 'w')
pwt.write("#!/bin/bash\n")
pwt.write("#SBATCH --job-name=preweighted\n")
pwt.write("#SBATCH --partition=sixhour\n")
pwt.write("#SBATCH --output=preweighted_%A.out\n")
pwt.write("#SBATCH --nodes=1\n")
pwt.write("#SBATCH --ntasks-per-node=2\n")
pwt.write("#SBATCH --constraint=intel\n")
pwt.write("#SBATCH --mem=30G\n")
pwt.write("#SBATCH --time=06:00:00\n")
pwt.write("module load Dir_Calc_Fluct\n")

pwt.write("cd $SLURM_SUBMIT_DIR\n")

pwt.write("rm LJD_init.out LJAnew_init.out LJAold_init.out\n")
pwt.write("rm CD_init.out CAnew_init.out CAold_init.out\n")
pwt.write("for i in {%d..%d..%d}; do cat FILES/$i/LJAold_init.out >> LJAold_init.out; done\n" % (inputparam.start_config, inputparam.end_config, inputparam.sep_config))
pwt.write("for i in {%d..%d..%d}; do cat FILES/$i/LJAnew_init.out >> LJAnew_init.out; done\n" % (inputparam.start_config, inputparam.end_config, inputparam.sep_config))
pwt.write("for i in {%d..%d..%d}; do cat FILES/$i/LJD_init.out >> LJD_init.out; done\n" % (inputparam.start_config, inputparam.end_config, inputparam.sep_config))
pwt.write("for i in {%d..%d..%d}; do cat FILES/$i/CAold_init.out >> CAold_init.out; done\n" % (inputparam.start_config, inputparam.end_config, inputparam.sep_config))
pwt.write("for i in {%d..%d..%d}; do cat FILES/$i/CAnew_init.out >> CAnew_init.out; done\n" % (inputparam.start_config, inputparam.end_config, inputparam.sep_config))
pwt.write("for i in {%d..%d..%d}; do cat FILES/$i/CD_init.out >> CD_init.out; done\n" % (inputparam.start_config, inputparam.end_config, inputparam.sep_config))
pwt.write("combine_weighted.py -etype flucts.special.inp -corr_name crp -dcorr_name dcrp -mol_name water\n")
pwt.write("combine_weighted.py -etype flucts.special.inp -corr_name frame -dcorr_name dframe -mol_name water\n")
pwt.write("touch .flag_comboflucts\n")
pwt.close()
