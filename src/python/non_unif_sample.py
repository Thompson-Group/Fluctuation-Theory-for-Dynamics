import numpy as np
from src.python.read_input import user_input

inputparam = user_input('input_file')

pstofs=1000.0


end=inputparam.nve_length+50.0
step = 0.0
f=open('time.dat','w')
g=open('real_time.dat','w')
while step < end:
    if step < 500:
        step = step + .01*pstofs
    else:
        step = step + inputparam.timestep*pstofs
    f.write("%d\n" % step)
    g.write("%s\n" % str(step/inputparam.sep_config))

