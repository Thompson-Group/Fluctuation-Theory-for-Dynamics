# Python Codes

This directory stores all the python codes that are used by the program to calculate various things. Here, I will briefly summarize everything that is in this directory, and which things are called when.

## Codes Directly Called by backbone.py 

• read\_input.py: This program reads the input file as a python class and it passes the information to the various python scripts.

• gen\_input.py: This program creates the above input file that is read by the above code. This is called automatically by backbone.py

• gen\_sub\_scripts.py: This program is called by the backbone.py script to generate all the submission scripts that are used by the code to submit all the jobs. In particular, this creates the job array that launches the jobs, it creates the segment creation array, it creates the segment combination code, and it creates the individual nve submission scripts.

• checkup.py: This program checks each of the directories in the FILES subdirectory and looks for the two files (log.lammps, corr\_complete) and if they don't exist it writes the file necessary to submit the individual nve run to recalculate it because it failed the first time. This is then automatically submitted if necessary.

• unif\_sample.py: This program writes out the correct times as real\_time.dat which is then used by the analysis codes to add them to the input file.

• non\_unif\_sample.py: This program writes out the correct times but now with non-uniform sampling based on a set of rules that are based on the information included in the nve file. This hasn't been incredibly successful, and should probably be avoided.

## Codes Used in the NVE Trajectories

• set\_msd\_calcs.py: This code creates a specific input file for the nve trajectory that tells each of the sub codes (including the fortran codes) what directory and files to read, as well as the simulation volume (for NPT simulations).


• grab\_press.py: This code looks at the lammps log file and grabs the pressures from the log file. Note: currently it assumes that the pressures are in the lammps log file as usecols=(13,14,15,16,17,18)) so taking up columns 14-19. This will be updated in a future release to actually seek the correct pressures.

• jump\_rot.py: This code calculates a large number of correlation functions (for molecules named "water" only). This calculates the hydrogen bond jump and frame correlation functions (calculate for the OO and OH frame - as well as C1-C3). It also calculates the hbond survival probability correlation function.

• cn\_rot\_calc\_special.py: This code calculates the frame reorientation for a single pair of atoms. This is useful if you have a single species that you want to know the reorientation of, OR, if you add non-molecular restraints to your system. Not called automatically, you must add this by hand to the run\_array.sh prior to the nve run.

• vel\_reselect.py: Code that is used for ion pairing simulations. This reselects the velocities of the ion pair by choosing them from the maxwell boltzmann distribution. Alpha code, use with care.

## Codes Called by Submission Scripts

• grab\_flucts.py: This reads the energies in from the lammps log files and outputs them to files based on the components defined in the flucts.inp file and puts them into the X\_init.out file.

• parse\_and\_combine.py: This reads the binary files outputted in the trajectory and writes new binary files that combine the first time that the code is run, replacing init\_segments.py. The second time it is run it calculates the derivatives and outputs block averaged results to the OUT/ folder, replacing combine\_segments.py. This last step can be rerun with modifications to the input\_file to change things like number of blocks, etc. 

• init\_segments.py (Obselete): This code reads in all the files in the different directories and then weights them by the fluctuations in energy by subdividing the problem into a user defined set of subproblems.


• combine\_segments.py (Obselete): This code reads in the subproblems from the previous code and calculates the derivatives and outputs the final derivative files.

• combine\_weighted.py: This is an alpha code that is used for deeper decompositions. By default this is turned off, but can be added in again by adding the flag .flag\_comboflucts to the directory.

## Codes  for Post-Analysis that are NOT Called Automatically

• msd\_fit.py: Calls the fitting of the diffusion coefficients and their derivatives.

• reor\_fit.py: Calls the fitting of the reorientation times and their derivatives. 

• jump\_fit.py: Fits the jump, frame, and integrates the theta contributions to the activation energy. Outputs in an easy to read manner. Currently an Alpha Code.

## Unused Codes:

calc\_cv\_cp.py and calc\_init\_energy.py are currently in development and unused.


