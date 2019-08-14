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
# reruns = [540, 543, 545, 547, 548, 550, 555, 560, 563, 564, 565, 568, 570, 572, 574, 575, 577,  579, 581, 582, 583, 586, 588, 592, 593, 595, 596, 598, 599, 601, 605,
#           606, 607, 608, 609-615, 617, 619, 621, 622, 623, 631, 636, 637, 639, 640, 642, 644, 645, 648, 649, 651, 653, 655, 657, 659, 661, 664, 665, 666, 668, 671,
#           672, 673, 674, 676, 677, 678,  679, 680, 682, 683, 684, 685, 687, 688, 689, 690, 691, 693, 697, 701, 704, 707, 708, 711, 712, 714, 716, 717, 721, 723, 725,
#           726, 727, 728, 729, 730, 731, 733, 735, 737,  741, 743, 745, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 761, 762, 764, 767, 769, 771, 775,
#           777, 780, 782, 783, 786, 788, 790]
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
tlist = [item for sublist in tlist for item in sublist]

# dd = [32, 39, 40, 41, 42, 43, 44, 45, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
#       42, 43, 44, 45, 46, 47, 48, 49]
# tlist = [22, 22, 22, 22,22,22,22,22, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
#          30,  30, 30, 30, 30, 30, 30, 30]

if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
    taskId = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    jobId = int(os.getenv('SLURM_JOB_ID'))
    print(taskId)
    # i = 10
else:
    # i = 12

s1 = str(dd[i-1])
t = str(tlist[i-1])

# open .m file
if "SLURM_JOB_ID" in os.environ:
    filepath = 'integrate_posterior.m'
    copyfile(filepath,'integrate_posterior_' + str(jobId) + '_' + str(taskId) + '.m')
    filepath = 'integrate_posterior_' + str(jobId) + '_' + str(taskId) + '.m'

else:
    filepath = '/Users/sarahfletcher/Documents/MATLAB/Repository_SDP/integrate_posterior.m'

with open(filepath, 'r') as file:
    data = file.readlines()

data[33] = 's1 = ' + s1 + ';' + '\n'
data[34] = 't = ' + t + ';' + '\n'

with open(filepath, 'w') as file:
    file.writelines( data )


with open(filepath, 'r') as file:
    data = file.readlines()
s1 = data[33][5:-2]
t = data[34][4:-2]
print(s1)
print(t)




