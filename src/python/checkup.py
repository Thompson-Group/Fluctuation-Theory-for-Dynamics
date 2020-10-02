#checkup.py This program just checks the log file of every sub simulation and makes sure that it exists - if it doesnt then it outputs the commands necessary to rerun the file.

import numpy as np
import os.path

with open("file_names") as f:
    for file in f:
        filepath=str("FILES/")+str(file).strip()+str("/log.lammps")
        filepath2=str("FILES/")+str(file).strip()+str("corr_complete")
        if os.path.isfile(filepath) and os.path.isfile(filepath2):
            output=1
        else:
            cd=str("cd FILES/")+str(file).strip()
            output=str("msub ")+str("nve.sh")
            back="cd ../../"
            print cd
            print output
            print back
        #DONE
    #DONE
#DONE
