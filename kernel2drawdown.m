function [drawdown] = kernel2drawdown(Q, kernel, t )
% Given previously calculatd kernel functions, pumping rates, and time
% period of interest, this funciton calculates the drawdown in each
% observation well using linear superposition.

% t is a vector which provides the time period over which the drawdown is
% calculated. If t = 4, that means we are only calculating the change in
% head between t = 3 and t = 4. If t = 1:4, we are calculating the entire
% drawdown since t = 1 to t = 4.

% The kernel functions provide MUST be for the corresponding vales of t.
% For example, if t = 3:5 are given, the kernel matrix must have 3 columns
% corresponding to time periods 3,4, and 5.

% Q is a [numTime X numPumpWells] matrix
% kernel is a [numObserve x length(t) x numParameterValues/runs x numPumpWells] matrix
% Drawdown is a [numObserve x length(t) x runs] matrix

% Check that Q and kernel have compatible dimensions
[Qtime, Qpump] = size(Q);
[~,ktime,~, kpump] = size(kernel);

if ~(Qpump == kpump)
    error('Different number pumping wells in Q and kernel')
elseif ~(ktime == Qtime)
    error('Different num time steps in Q and kernel')
end

[ddsize1, ddsize2, ddsize3,numPumpWells] = size(kernel);
drawdownKernelComponent = zeros([ddsize1, ddsize2, ddsize3]);
for i = 1:numPumpWells
    for j = 1:length(t)
        drawdownKernelComponent(:,j,:) = drawdownKernelComponent(:,j,:) + kernel(:,j,:,i) * Q(j,i);
    end
end
drawdown = cumsum(drawdownKernelComponent,2);
