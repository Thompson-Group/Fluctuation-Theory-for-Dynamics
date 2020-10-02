This code is copyright October 2020 by Ezekiel Piskulich and Ward Thompson, University of Kansas. Users are free to share/modify this code as desired; however, this code shall not be incorporated into proprietary software without prior permission from the author of this code. All rights reserved. Any questions about this repository can be directed at zeke.piskulich at gmail dot com. 

# Fluctuation Theory for Dynamics (DCM)

This code is a general use code developed to calcualte derivatives (T,P) of time correlation functions (TCFs) from Molecular Dynamics simulations. The code thus far includes support for:

    - Mean-Squared Displacement

    - Reorientation (C1, C2)

    - Jump Reorientation 

    - Viscosity (Shear)

    - Reactive Flux (Flux-Side)

Adding new correlation functions is relatively simple, you just need to add a fortran code and then the subsequent calls to the submission script generation python code (gen\_sub\_scripts.py).

## Installation Instructions
Upon first installation it is necessary to make the program, which is a two step process. 

1) Step One: Make Program
    - Here you need to set information about your machine and your compiler as well as any compilation flags you want to include.
    - Edit the makefile in the main directory.
    -Then type:
    '''
    make
    '''
    - By default this adds a modulefile for loading on a cluster that loads by the cluster and adds an appropriate line to your bash\_profile. If your cluster does not use modules, delete this line and instead add the bin/ folder to your PATH environment.

2) Step Two: Update Header File:
    - Whatever option you chose for MACHINE in your MAKEFILE needs to have a $MACHINE\_header.dat in src/dependencies. Example header files are included in this directory for the KU Community Cluster and the NERSC Cori Computer.

3) Step Three: If not running on the CRC at KU, minor code modifications need to occur.
    - In src/python/gen_sub_scripts.py need to change all SBATCH directives to appropriate ones for your cluster.
    - You must change src/shell/grabfluctssub.sh to match your system settings.
    - You must change the header of src/gen_sub_scripts.py to match your system settings. 
    - This code depends primarily on slurm commands, if you aren't using slurm, significant modifications may be necessary.

## How to Use the Program
On the simplest level, running the program is just repeatedly typing:
'''
backbone.py
'''
Note: backbone.py can be run with the option -dryrun 1 to avoid automatic submission of jobs.

This script runs each step of the code and checks for completion, and then spits out what you should do next.
This script is split up into steps that are each their own python function.
These steps are:
1) MACHINENAME: reads the machine name and the DCM path

2) INSTRUCT: Gives initial instructions and lets you generate input file

3) TRAJ: Checks that the long trajectory has been completed

4) NVE: Prepares the short trajectories

5) NVECOMPLETE: Makes sure you don't move on until trajectories are complete

6) CHECKUP: Checks that no runs failed prematurely

7) GRABFLUCTS: Grabs fluctuation info from the log files of each nve run.

8) SEGARRAY: Reads in all the binary files from the nve trajectories and combines them into one binary file.

9) COMBINESEG: Combines the segments into block averages and total avererages.


