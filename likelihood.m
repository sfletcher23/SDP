% NN info
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
gwParam.startingHead = 337.143;

t = 10;
logk = 0:.01:5;
logs = -12:.01:-11;
s1 = 10;
llh = zeros(length(logk), length(logs), length(s1));
for i = 1:length(logk)
    for j = 1:length(logs)
        for k = 1:length(s1)
    
            
            s = log(1E-5);

            % Calculate likelihood using model
            input = [logk(i); logs(j); 365*t];
            drawdown_t_current = netscript(input, gwParam);
            y = drawdown_t_current;
            u = s1(k); 
            llh(i,j,k) = normpdf(y, u, L_sigma);
        end
    end
end

