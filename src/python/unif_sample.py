import numpy as np
from src.python.read_input import user_input

# Call input file class.
inputparam = user_input('input_file')

pstofs = 1000.0

step = 0.0
f=open('time.dat','w')
g=open('real_time.dat','w')
while step < inputparam.nve_length:
    f.write("%d\n" % step)
    g.write("%s %s\n" % (str(step/pstofs),1.0))
    step = step + inputparam.timestep*pstofs
