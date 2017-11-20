function [y1] = myNeuralNetworkFunction_54215(x1, gwParam)
%MYNEURALNETWORKFUNCTION_54215 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 16-Nov-2017 21:11:36.
% 
% [y1] = myNeuralNetworkFunction_54215(x1) takes these arguments:
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
b1 = [-0.47716324761905859475;-1.497468837221651583;2.5599702393927010746;0.48642067733346866509;0.45409106325491871425;-1.5628432113303234452;-3.767957545057211366;3.3705286872343194204];
IW1_1 = [2.0962518178963893689 0.29847969200409657509 -0.60022978282897665459;1.0348102855052538906 0.13255846468020637818 2.3545070226680318015;-0.6818820992090151023 -1.1533320115157581487 -0.36779333973885464104;1.5918239293761551334 1.4536270346034827483 1.7210424676260553944;0.83144594363537493997 2.3739996807513428223 0.0082386270134222790162;-2.1367478989420547464 -0.33580856899649169645 0.3828812271971077652;-3.3386109527155363175 -0.52714968647848536865 0.73172927039710888852;1.4563561831924245471 0.20387728986066377956 2.1961896832659384415];

% Layer 2
b2 = 2.6453650715535985682;
LW2_1 = [0.15355959848282341174 0.043377921641841718503 0.067682267618461106062 0.0087922134144255356669 -0.0073985493276802795359 -1.9465987216703057872 2.72724980208316925 -1.0775867608579230428];

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