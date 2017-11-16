function [y1] = myNeuralNetworkFunction_49010(x1, gwParam)
%MYNEURALNETWORKFUNCTION_49010 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 13-Nov-2017 11:43:19.
% 
% [y1] = myNeuralNetworkFunction_49010(x1) takes these arguments:
%   x = 3xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% time = x1(3,:);
% time = time - 365;
% indexZeroTime = time < 0; 
% x1(3,:) = time;

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-4.18947799393779;-12.2047587917482;3.65000009536743];
x1_step1.gain = [0.243621113364753;1.24386853941395;0.000182706957739474];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.9875103882231148766;-1.9590839918992666302;-0.64639608688157057514;-0.35126327062947304558;0.51044278075696603025;-1.9921670973538212479;2.6273507381805556449];
IW1_1 = [-0.238784209728525737 -0.079954815010411517395 -2.6602475277663049091;-0.21938961412831622932 -0.59785653185646470131 -2.3822186525520634426;1.0287201608762936988 2.2701557414405884039 -0.63431209993340698094;2.0003629858216540605 0.29006341957849257618 -0.17647174955011235564;1.5800840842876140435 -0.032181896915685216654 1.7845990099869968315;-0.77034602277128416681 0.40777391474736518484 -2.3644848662261206407;3.2654177120727139183 0.5135558880220933986 -0.42899950227576849615];

% Layer 2
b2 = [1.8864116283527117002;1.497338837532039113;0.32787358983506686982;-0.72232786297952100707;0.39401959207860670631;-0.9866480307799784244;-2.050394708071941352];
LW2_1 = [-0.81841953142443246705 0.88263578413856169647 -0.43510953480525849102 -0.23546067370622073001 0.38041788051410169169 0.22672742527038258764 1.2290660576706355922;-1.2708750714789915914 0.39828602679713348556 0.30704643984483059871 0.60313729448205555261 0.011429170169056833417 -0.67325422244851496067 -0.78886579454007210721;0.069254236119569684282 0.22888515814428753159 0.29228343732858008686 -1.4258829686370972922 -0.46038112642858236567 0.074479613202645719716 1.7717777952978699396;-0.88966817996479852226 -0.82981083271654976219 -0.22115210928982725891 0.56897813417453890761 -0.97905910346695168567 -0.47976583741737455213 1.6277674963576549416;0.82160883787154237456 -0.16550196296265148876 -0.067482571988763714876 0.94716714891523801079 0.30408040686113763096 0.12810465844390861445 1.4167904379836795403;-0.48015227316366471122 -0.24143391367285660376 -0.03646313695344771677 -0.5994062865812647134 0.28871981517586525756 -0.25972977751307835259 -0.94536111700346037168;-0.75628171570429403303 0.17701118354978082059 -0.85316566049236575164 0.72237184081229355304 0.33950651810818605059 -0.20832650237408961713 -1.4350021052793748844];

% Layer 3
b3 = [1.9403539838075347657;-1.1075068297339349943;0.54578720974180161551;0.54030689486818017642;-1.0232093906726467125;1.0347749477175713917;-1.6799794530916722923];
LW3_2 = [-0.3224403002930181783 0.77174824304699596311 -1.046202988182223903 1.0658425399591557881 0.058023249567605281185 0.41731297195474781336 0.44731237938321150027;0.57342990933251414276 0.69249230185981003327 1.7461212977382514211 1.8498866086717338231 -0.81774737928760399264 -1.1333658999422475677 -1.0906523194976216473;-0.82168668045826298574 0.58571717949704726891 -0.28460251644695039319 -0.15017629618588582141 1.007208799608054095 -1.2156957545097835105 0.35405471033391588964;1.2162078457133389531 1.2242617960905308738 0.29903581053669886503 -0.28565396559098477081 1.1726153877146212601 -0.74947459264319993277 -0.45024193780718846369;-1.2148642262425262128 -0.90449163673511867589 0.46776081637927929835 -0.56800090075846354498 0.011762229003052722051 -0.36870972103440818879 1.1755644700683960391;0.27221885756754388153 0.65073480617098200529 0.31248498706500393851 1.3769617883452096585 -0.81888703842361920682 1.2722244282281076622 0.53326634833138908398;0.14505814565099089952 -0.49396055596370835872 0.1708268540364197563 -1.1236397378632732025 -1.1079511168355495876 0.85737501065401999778 0.58649527863881456557];

% Layer 4
b4 = 0.809633099837997805;
LW4_3 = [0.072579264157358588938 -0.99386699273720058301 0.83074627907777587321 -0.084196372405547739715 -0.54868139260163750759 -1.0156623187454756341 -0.60470792207729351375];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00215538909597944;
y1_step1.xoffset = -587.998901367188;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);

% Layer 3
a3 = tansig_apply(repmat(b3,1,Q) + LW3_2*a2);

% Layer 4
a4 = repmat(b4,1,Q) + LW4_3*a3;

% Output 1
y1 = mapminmax_reverse(a4,y1_step1);


% Convert to drawdown
y1 = gwParam.startingHead - y1;

% Set initial head to starting
%y1(indexZeroTime) = 0;

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