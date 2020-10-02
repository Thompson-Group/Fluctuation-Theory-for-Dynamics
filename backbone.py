#!/usr/bin/env python

import numpy as np
import os,sys,subprocess
import src.python.path

def step_machinename():
    """
    This step checks if the machine.name file exists in the directory in which you are running.
    This file specifies the machine that you are running on, in a single line.
    This machine should match the information you provide in the makefile.

    This is set by the import statement up above, and works only if the code has been made using the supplied makefile.
    """
    exists = os.path.isfile('machine.name')
    if exists:
        f = open('machine.name', 'r')
        mach_name = f.readline().strip()
        homepath = f.readline().strip()
        f.close()
    else:
        print("->machine.name is not set!")
        sys.exit()
    return mach_name, homepath

def step_instruct():
    """
    This step gives the basic instructions for completing a run. It makes sure you provide the 
    two needed files, input_file, and flucts.inp
    It also lets you call the input generation code.
    """
    flag='.flag_instruct'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Instructions Flag exists.")
    else:
        print("Instructions Flag doesn't exist")
        print("->Please create input_file")
        print("->Please create flucts.inp")
        if dryrun == 0:
            print("->Alternatively, you can run the input_generator to create input_file")
            inpresponse = str(input("Would you like to do so? y/n\n"))
            print("You have selected: %s" % inpresponse)
            if 'y' in inpresponse:
                import src.python.gen_input
            elif 'Y' in inpresponse:
                import src.python.gen_input
        if not os.path.isfile("input_file"): sys.exit("Error: the file input_file missing")
        if not os.path.isfile("flucts.inp"): sys.exit("Error: the file flucts.inp missing")
        if not os.path.isfile("corr.funcs"): 
            print("corr.funcs not found: generating automatically. Feel free to add correlations")
            print("Currently has msd, cn as options.")
            print("Other options: fsc_f, fsc_b, shear, theta, frame, crp")
            corrs=["msd","cn"]
            with open("corr.funcs",'w') as f:
                for item in corrs:
                    f.write("%s\n" % item)
        open(flag, 'a').close()
        print("Instructions step complete, please run backbone.py again.")
        sys.exit()
    return


def step_traj():
    """
    This is the code that checks that an NVT or NPT trajectory has been completed.

    This can work with a lammps run (log.lammps) or a cp2k run (cp2k.log).
    """
    flag='.flag_traj'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Trajectory Run Flag Exists.")
    else:
        print("Trajectory Run Flag Missing.")
        f = open('mol_names','w')
        for item in inputparam.molec:
            f.write('%s\n' % item)
        f.close()
        if inputparam.prog == "LAMMPS":
            log_exists = os.path.isfile('log.lammps')
            if log_exists:
                print("->Trajectory found")
                open(flag, 'a').close()
            else: 
                print("-->At this point, you should run the base trajectory")
                print("-->Example NVT and NPT input files are located in src/dependencies")
                sys.exit()
        elif inputparam.prog == "CP2K":
            log_exists = os.path.isfile('cp2k.log')
            if log_exists:
                print("->Trajectory found")
                if not os.path.exists("RESTART"):
                    os.makedirs("RESTART")
                subprocess.call("mv *.restart RESTART", shell=True)
                print("->What is the identifier for your restart files?")
                identifier=str(input("-->i.e. anything before the # in the .restart files\n"))
                if os.path.isfile('RESTART/'+identifier+str(inputparam.start_config)+".restart"):
                    print("->Renaming restart files")
                    for i in range(inputparam.start_config,inputparam.end_config,inputparam.sep_config):
                        os.rename('RESTART/'+identifier+str(i)+'.restart','RESTART/restart.'+str(i))
                    open(flag, 'a').close()
                else:
                    print("-->Improper identifier, try again.")
                    print("-->If you would like to skip, please type touch .flag_traj")
            else:
                print("->At this point, you should run the base trajectory")
                print("->Example NVT and NPT input files are located in src/dependencies")
                sys.exit()
        else:
            print("Incorrect program type")
            sys.exit()
    return 

def step_nve():
    """
    This step makes the file_names file, the FILES directory, and then calls a couple of python
    codes to create submission scripts etc
    """
    if not os.path.isfile('in.nve') and not os.path.isfile('in.nve.cp2k'):
        print("Need to put an in.nve or in.nve.cp2k file into directory")
        sys.exit()
    for i in range(inputparam.num_molecs):
        if inputparam.molec[i] != "water" and not os.path.isfile(str(inputparam.molec[i])+'.txt'):
            print("Need to create %s.txt" % inputparam.molec[i])
    flag='.flag_nve'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("NVE Flag Exists.")
    else:
        print("NVE Flag Missing.")
        print("->Setting up and running NVE trajectories")
        f = open("file_names", 'w')
        # Writes the filenames
        for i in range(inputparam.start_config, inputparam.end_config+inputparam.sep_config, inputparam.sep_config):
            f.write("%d\n" % i)
        f.close()
        if not os.path.exists("FILES"):
            os.makedirs("FILES")
        if corr_override == "None": import src.python.gen_sub_scripts
        else: import src.python.gen_sub_override
        if inputparam.timestep == 'FALSE':
            import src.python.non_unif_sample
        else:
            import src.python.unif_sample
        print("->NVE Trajectories are ready to be submitted as sbatch run_array.sh")
        open(flag, 'a').close()
        sys.exit()
    return

def step_nvecomplete():
    """
    This step exists to make sure that you don't go to the next real step 
    if simulations are still running. It basically just makes you double check
    that all jobs have finished.
    """
    flag='.flag_nvecomplete'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("NVE Completion Flag Exists.")
    else:
        print("NVE Completion Flag Missing.")
        print("->Note: you have to create this flag yourself")
        print("->Check the queue - make sure ALL jobs have finished")
        print("->Then, once you are satisfied type 'touch .flag_nvecomplete'")
        sys.exit()
    return

def step_checkup():
    flag='.flag_checkup'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("NVE Trajectory Checkup Flag Exists.")
    else:
        print("NVE Trajectory Checkup Flag Missing.")
        print("-> Running Checkup")
        f = open('rerun', 'w')
        corr = np.genfromtxt('corr.funcs',dtype=str)
        chk = corr[-1]
        if chk == "frame": chk = chk+"c3"
        elif chk == "cn": chk = "c2"
        elif chk == "msdpy": chk = "msd_py"
        for i in range(inputparam.start_config, inputparam.end_config+inputparam.sep_config, inputparam.sep_config):
            if (i%1000000) == 0:
                print("-->%d has been checked" % i)
            if inputparam.cab == "TRANSPORT" and not os.path.isfile("FILES/%d/log.lammps" % (i)):
                print("%d does not exist" % i)
                f.write("mkdir FILES/%d; cp in.nve FILES/%d; cp nve.sh FILES/%d; cd FILES/%d; sed -i -e 's@AAA@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i,i))
            elif inputparam.cab == "TRANSPORT" and not os.path.isfile("FILES/%d/%s_%d_%s.pckl" % (i,chk,i,inputparam.molec[0])):
                print("%d does not exist" % i)
                if inputparam.prog == "LAMMPS":
                    f.write("mkdir FILES/%d; cp in.nve FILES/%d; cp nve.sh FILES/%d; cd FILES/%d; sed -i -e 's@AAA@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i,i))
                elif inputparam.prog == "CP2K":
                    f.write("mkdir FILES/%d; cp in.nve.cp2k FILES/%d; cp nve.sh FILES/%d;  cd FILES/%d; sed -i -e 's@AAA@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i,i))
            elif inputparam.cab == "IONPAIRING" and not os.path.isfile("FILES/%d/fsc_f_%d_%s.dat" % (i,i,inputparam.molec[0])):
                print("%d does not exist" % i)
                if inputparam.prog == "LAMMPS":
                    f.write("mkdir FILES/%d; cp in.nve FILES/%d; cp nve.sh FILES/%d; cd FILES/%d; sed -i -e 's@AAA@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i,i))
                elif inputparam.prog == "CP2K":
                    f.write("mkdir FILES/%d; cp in.nve.cp2k FILES/%d; cp nve.sh FILES/%d;  cd FILES/%d; sed -i -e 's@AAA@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i,i))
        f.close()
        if os.path.isfile('rerun') and os.path.getsize('rerun'):
            f = open('rerun','r')
            count = 0
            for line in f:
                count += 1
            f.close()
            print("-->There are %d missed jobs" % count)
            if count < 100:
                subprocess.call("bash rerun", shell=True)
                print("-->Submitting %d jobs" % count)
            else:
                print("-->Too many jobs, please check and submit by bash rerun")

        print("->Wait until remaining NVE trajectories run")
        print("->If this step doesn't work the first time")
        print("->Check your input files")
        open(flag, 'a').close()
        sys.exit()
    return

def step_comboflucts():
    """
    This step calls the enhanced decompositions
    Eventually this flag will become depreciated.
    By default, this step is skipped
    """
    flag='.flag_comboflucts'
    flag_exists = os.path.isfile(flag)
    sp_exists = os.path.isfile("flucts.special.inp")
    if not flag_exists:
        print("Combo Flucts Flag Exists,skipping.")
    else:
        print("Combo Flucts Flag Missing")
        if "water" not in inputparam.molec:
            open(flag, 'a').close()
        else:
            if not sp_exists:
                print("flucts.special.inp doesn't exist")
                sys.exit()
            else:
                os.system('sbatch comboflucts.sh')
                print("-> Wait until comboflucts has finished before running the next step!")
                print("-> This step will write the appropriate flag upon completion.")
                sys.exit()
    return


def step_grabflucts():
    """
    This step calculates the energy fluctuations in the TOTAL energies by looking at each directory
    and combining them.
    """
    flag='.flag_grabflucts'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Grab Flucts Flag Exists.")
    else:
        print("Grab Flucts Flag Missing.")
        if not os.path.isfile("flucts.inp"):
            print("flucts.inp missing")
            print("Please generate.")

        #subprocess.call("cp %s/src/shell/grabfluctsub.sh ./" % homepath, shell=True)
        if dryrun == 0: os.system("sbatch grabfluctsub.sh")
        print("->Grab flucts is running as a job")
        print("->Please be patient - and check that grab_flucts.py is completed before continuing")
        print("->Grab Flucts has Finished!")
        open(flag, 'a').close()
        sys.exit()
    return

def step_segarray():
    """
    This divides the trajectories into more manageble chunks that speed up things
    like block average calculations and puts them into the SEG directory
    """
    flag='.flag_segarray'
    flag_exists = os.path.isfile(flag)
    # Creates electrostatics
    if os.path.isfile('ewald_init.out') and os.path.isfile('coul_init.out'):
        ewald,coul = np.genfromtxt('ewald_init.out',usecols=0),np.genfromtxt('coul_init.out',usecols=0)
        elec = np.add(ewald,coul)
        np.savetxt('elec_init.out',np.c_[elec])
        with open('flucts.inp','r') as f:
            lines = f.readlines()
            segflag=0
            for line in lines:
                if "elec" in line:
                    segflag = 1
        if segflag != 1:
            f = open("flucts.inp",'a')
            f.write("elec 0\n")
            f.close()
    if flag_exists:
        print("Segment Creation Array Flag Exists.")
    else:
        print("Segment Creation Array Flag Missing.")
        print("->Running Segement creation calc.")
        if not os.path.exists("SEG"):
            os.makedirs("SEG")
        if dryrun == 0: os.system("sbatch init_segments.sh")
        print("->Segment Creation Array Submitted")
        print("->Please wait until all jobs are complete")
        open(flag, 'a').close()
        sys.exit()
    return

def step_combineseg():
    """
    This combines the segments from the previous step and block averages
    etc to come up with final answers
    """
    flag='.flag_combineseg'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Segment Combination Flag Exists.")
    else:
        if not os.path.exists("OUT"):
            os.makedirs("OUT")
        print("Segment Combination Flag Missing.")
        print("->Running segment combination")
        if dryrun == 0: os.system("sbatch combine_segments.sh")
        print("->Segment Combination has finished!")
        open(flag, 'a').close()
        sys.exit()
    return


if __name__ == "__main__":
    mach_name, homepath = step_machinename()

    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('-dryrun', default=0,type=int,help="Dry run flag - [0] Submit Scripts as Jobs, [1] Creates, but doesn't submit")
    parser.add_argument('-corr_override', default="None", type = str, help="Flag that specifies correlation function commands to include in calculation files, if necessary")
    args = parser.parse_args()

    dryrun = args.dryrun
    corr_override = args.corr_override
    if corr_override != "None":
        f = open(".flag_override", "w")
        f.write(corr_override)
        f.close()


    step_instruct()
    from src.python.read_input import user_input
    inputparam = user_input("input_file")
    step_traj()
    if not os.path.exists("logs"):
        os.makedirs("logs")
    step_nve()
    step_nvecomplete()
    step_checkup()
    step_grabflucts()
    step_comboflucts()
    step_segarray()
    step_combineseg()

