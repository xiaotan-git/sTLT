function [grid,data_new] = build_unicycleCBF_helper_data(grid,data_origin,backward_time,vRange,wMax)
%BUILD_UNICYCLECBF_HELPER Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    vRange = [-1 1];
    wMax = 1;
end

t0 = 0; Nt = 41;
tMax = backward_time; % todo: what if tset_low = tb_low
tau = linspace(t0,tMax,Nt);

%% reachability problem parameters
uMode = 'min';
dUni = Plane([0,0,0],wMax,vRange);

% Put grid and dynamic systems into schemeData
schemeData.grid = grid;
schemeData.dynSys = dUni;
schemeData.accuracy = 'high'; %set accuracy
schemeData.uMode = uMode;

%% Compute value function

%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 1;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.quiet = true;

% uncomment if you want to see a 2D slice
%HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
%HJIextraArgs.visualize.plotData.projpt = [0]; %project at theta = 0
%HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%[data, tau, extraOuts] = ...
% HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
[data, ~, ~] = ...
  HJIPDE_solve(-data_origin, tau, schemeData, 'minVOverTime', HJIextraArgs); % careful with the reverse sign

% now we flip the time and change the sign back
dataCBF = -flip(data,4);


data_new = squeeze(dataCBF(:,:,:,1));
end

