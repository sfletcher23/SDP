import os
import numpy as np
from shutil import copyfile


time = np.arange(1, 30, 1)
drawdown_max = [ 71.8024349691912,	98.9996385259814,	121.941450801024,	141.343687632105,	157.856744357095,	172.048031095633,	184.398168603787,
                         195.305381025659,	205.094059319596,	214.024853008525,	222.304709021687,	230.096011944031,	237.524453528662,	244.685548197193,
                         251.649884766433,	258.467318232683,	265.170400521298,	271.777451521216,	278.295783299501,	284.725670246833,	291.065599430146,
                         297.318968271482,	303.501574870096,	309.648071618181,	315.814627673233,	322.075373217118,	328.512369315768,	335.202000802936,
                         342.202671675442,	349.547910713690]
drawdown_min = [7.90152696173453,	10.8863953854826,	13.3341292794019,	15.3663710702284,	17.0783085592299,	18.5440767178084,	19.8211011075899,
                        20.9535357329021,	21.9749402601961,	22.9103236238453,	23.7776598095863,	24.5889611257012,	25.3509773897663,	26.0655789803655,
                        26.7298813924321,	27.3361843409086,	27.8718371562427,	28.3192122981390,	28.6560732294635,	28.8567461499751,	28.8945904598435,
                        28.7461885505084,	28.3972642466096,	27.8494647933049,	27.1259694602479,	26.2731049674530,	25.3557618531984,	24.4468338880858,
                        23.6140170959438,	22.9088843642109]
dd = []
tlist = []
numRums = 0
for i in range(len(time)):
    ceil = min(np.ceil(drawdown_max[i]), 50)
    statesThisTime = np.arange(np.floor(drawdown_min[i]), ceil, 1)
    numStatesNow = np.size(statesThisTime)
    numRums = numRums + numStatesNow
    dd.append(statesThisTime)
    tlist.append([i+1] * numStatesNow)
dd = [item for sublist in dd for item in sublist]
tlist = [item for sublist in t for item in sublist]

if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
    i = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    print(i)
else:
    i = 0

s1 = str(dd[i])
t = str(tlist[i])


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




