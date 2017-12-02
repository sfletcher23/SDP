% Calculate head in t+1 period from samples

taskID = getenv('SLURM_ARRAY_TASK_ID');
jobID = getenv('SLURM_JOB_ID');

filename = strcat('samples_', num2str(jobID), '_',  num2str(taskID), '.mat');
load(filename)

% NN info
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
gwParam.startingHead = 337.143;

% Calculate head at next period for each sample
input = [sample_logk; sample_logs; repmat(365*t, size(sample_logk))];
drawdown_current = netscript(input, gwParam);
input = [sample_logk; sample_logs; repmat(365*(t+1), size(sample_logk))];
drawdown_next = netscript(input, gwParam);
dd = drawdown_next - drawdown_current;

save(filename, 'drawdown_current', 'drawdown_next', 'dd')
