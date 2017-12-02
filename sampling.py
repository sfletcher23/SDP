import scipy.stats as st
import numpy as np
import scipy.io as io
import scipy.interpolate as intp
import scipy.integrate as integrate
import matplotlib as ml
import matplotlib.pyplot as plt
import os

num_sample = 10
plotOn = True
if "SLURM_JOB_ID" in os.environ:
    plotOn = False

# open .m file
if "SLURM_JOB_ID" in os.environ:
    filepath = 'integrate_posterior.m'
else:
    filepath = '/Users/sarahfletcher/Documents/MATLAB/Repository_SDP/integrate_posterior.m'

with open(filepath, 'r') as file:
    data = file.readlines()
s1 = data[30][5:-2]
t = data[31][4:-2]
print(s1)
print(t)

# Load estimated pdf values
data = io.loadmat('posterior_samples.mat')
p = data['norm_p']

S_lower = 6.09e-6
S_upper = 2.2e-5
K_lower = 0.0001
K_upper = 25

# Define function that interpolates in order to get an arbitrary pdf value for joint distrbution
k = np.arange(np.log(K_lower),np.log(K_upper), 0.01)
s = np.arange(np.log(S_lower),np.log(S_upper), 0.01)
joint_ks = intp.interp2d(k,s,np.transpose(p), kind='cubic')
# total_p = integrate.dblquad(joint_ks, s[0], s[-1], lambda x: k[924], lambda x: k[-1], epsabs=1.49e-05, epsrel=1.49e-05)
# print(total_p[1])


# Integrate to get marginal distribution for S
p_k = np.zeros(len(s))
for i in range(len(s)):
    print(i)
    [p_k[i], er] = integrate.quad(joint_ks,np.log(K_lower), np.log(K_upper), args=s[i], epsabs=1.49e-06, epsrel=1.49e-06, limit=200 )
marg_s = intp.interp1d(s, p_k, kind='cubic')
# Check that marginal integrates to 1
total_p = integrate.quad(marg_s, s[0], s[-1], limit=200)
print(total_p[0])
if abs(total_p[0] - 1) > 0.01:
    error('Marginal dist for S not valid')
else:
    np.save('sample_data', marg_s, joint_ks)
# plot marginal
if plotOn:
    plt.plot(s,marg_s(s))
    plt.show()


# Calculate conditional of K on S
p_k_s = np.zeros([len(k), len(s)])
for i in range(len(k)):
    for j in range(len(s)):
        p_k_s[i,j] = joint_ks(k[i], s[j]) / marg_s(s[j])
cond_k_s = intp.interp2d(k,s, np.transpose(p_k_s))
np.save('sample_data', marg_s, joint_ks, cond_k_s)


# Sample from marignal S
class my_marg_s(st.rv_continuous):
    def _pdf(self, s):
        return marg_s(s)
marg_s_dist = my_marg_s(name='marg')
marg_s_dist.a = s[0]
marg_s_dist.b = s[-1]
r = np.random.rand(num_sample)
sample_s = marg_s_dist.ppf(r)

# For each marginal S sample, sample K from conditional
sample_k = np.zeros(num_sample)
for i in range(len(sample_s)):
    class my_cond_k_s(st.rv_continuous):
        def _pdf(self,k):
            return cond_k_s(k,sample_s[i])

    cond_ks_dist = my_cond_k_s(name='cond')
    cond_ks_dist.a = k[0]
    cond_ks_dist.b = k[-1]
    r2 = np.random.rand(1)
    sample_k[i] = cond_ks_dist.ppf(r2)

jobId = os.environ.get('SLURM_JOB_ID', [])
outputDic = dict(zip(['sample_logs', 'sample_logk', 's', 't'], [sample_s, sample_k, s1, t]))

if "SLURM_JOB_ID" in os.environ:
    jobId = os.getenv('SLURM_JOB_ID')
    if "SLURM_ARRAY_TASK_ID" in os.environ:
        taskID = os.getenv('SLURM_ARRAY_TASK_ID')
        filename = 'samples_' + str(jobId) + '_' + str(taskID)
    else:
        filename = 'samples_' + str(jobId) + '_' + 'dd' + str(s1) + '_t' + str(t)
else:
    filename = 'samples_' + 'dd' + str(s1) + '_t' + str(t)

io.savemat(filename, outputDic)

