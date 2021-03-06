function [y1] = myNeuralNetworkFunction_54212(x1, gwParam)
%MYNEURALNETWORKFUNCTION_54212 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 16-Nov-2017 23:19:41.
% 
% [y1] = myNeuralNetworkFunction_54212(x1) takes these arguments:
%   x = 3xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-4.24284330816697;-12.2060149616686;3.65000009536743];
x1_step1.gain = [0.237653620258443;1.24278963206423;0.000182706957739474];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.0421746479754498971;0.66603592075575079878;0.74060755872876593564;-2.63184078110013564;3.05505765893768233;-2.2171377779171352884;-2.2548167523574029758];
IW1_1 = [1.1631581247706452853 -0.55618332766519362753 2.294806849779644331;1.4380466353368577259 0.22284383701069981898 -0.23906067227417321375;-1.4469289129248896764 2.1263756651479126525 -0.43516077674093928129;-2.2785735944906160277 -0.35756638695736780997 0.54325900361173062869;1.3561074986162520162 0.19242101150414064881 1.9253882149498495746;-1.4232894255908463688 0.89163557733666087834 1.7130372333070726931;-3.4638784864566156685 -0.53989111553671542687 0.65065379831729208693];

% Layer 2
b2 = 2.9072312454492146294;
LW2_1 = [0.024002222327702060073 0.98743688249001415613 -0.0057058973347592269881 2.7424932773887578108 -1.1432805353973614082 0.013572723175363102591 -1.0199311203944192439];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0021516005056734;
y1_step1.xoffset = -587.99462890625;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);

% Convert to drawdown
y1 = gwParam.startingHead - y1;
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
