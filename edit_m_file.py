import os
import numpy as np
from shutil import copyfile

sgw = np.arange(20, 30, 1)
time = np.arange(5, 10, 1)

if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
    i = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    print(i)
else:
    i = 0

s1 = str(sgw[i])
t = str(time[i])


# open .m file
if "SLURM_JOB_ID" in os.environ:
    filepath = 'integrate_posterior.m'
    copyfile(filepath,'integrate_posterior_' + str(i) + '.m')
    filepath = 'integrate_posterior_' + str(i) + '.m'

else:
    filepath = '/Users/sarahfletcher/Documents/MATLAB/Repository_SDP/integrate_posterior.m'

with open(filepath, 'r') as file:
    data = file.readlines()

data[30] = 's1 = ' + s1 + ';' + '\n'
data[31] = 't = ' + t + ';' + '\n'

with open(filepath, 'w') as file:
    file.writelines( data )


with open(filepath, 'r') as file:
    data = file.readlines()
s1 = data[30][5:-2]
t = data[31][4:-2]
print(s1)
print(t)




