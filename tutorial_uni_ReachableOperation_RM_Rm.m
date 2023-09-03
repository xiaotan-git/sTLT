% clc;
% unicycle setting
vRange = [-1 1];
wMax = 1;
timeInterval = [1,1.5];

X1.c = [-2;2]; X1.r = 1; 
[grid,data_original] = build_unicycleCBF_helper_cylinder([-4;-4;-pi],[4;4;pi],X1.c,X1.r);
Nt = 100;

grid_min = grid.min; grid_max = grid.max; N = grid.N;

pdDims = 3;
% g_with_time = createGrid([grid_min; tb_low], [grid_max;tset_low], [N; Nt], pdDims);

t0 = timeInterval(1);
tMax = timeInterval(2); 
tau = linspace(t0,tMax,Nt);

%% reachability problem parameters
dUni = Plane([0,0,0],wMax,vRange);

% Put grid and dynamic systems into schemeData
schemeData.grid = grid;
schemeData.dynSys = dUni;
schemeData.accuracy = 'high'; %set accuracy


%% Compute value function

%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 1;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.quiet = false;

% uncomment if you want to see a 2D slice
%HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
%HJIextraArgs.visualize.plotData.projpt = [0]; %project at theta = 0
%HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%[data, tau, extraOuts] = ...
% HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
% [data, ~, ~] = ...
%   HJIPDE_solve(data0, tau, schemeData, 'minVOverTime', HJIextraArgs);

% Define the two reachable problems
% RM(S,[a,b]), which corresponds to F[a,b]S
% to solve RM(S,[a,b]) we need to minimize over t and minimize over u
uMode = 'min'; data0  = -data_original; a = 2; b =3;
tau_F = linspace(0,b,Nt);
schemeData.uMode = uMode;
[data_F, ~, ~] = ...
  HJIPDE_solve(data0, tau_F, schemeData, 'zero', HJIextraArgs);
data_RM = -squeeze(data_F(:,:,:,end)); %tau

figure(2)
    clf
    h0 = visSetIm(grid, data_F(:,:,:,1));
    h0.FaceAlpha = .3;    hold on
    h = visSetIm(grid, data_F(:,:,:,end));
    h.FaceAlpha = .3;
    title('The reachable set at the end');
    x0 = [0;0;0];
    if eval_u(grid,data_RM,x0)>0
        plot3(0,0,0,'*r');
    else
        plot3(0,0,0,'g');
    end
    hold off
    


% Rm(S,[a,b]): all states from which no matter what input signal u ∈ U is applied, 
% the system can reach the target set S at some time instant in [a,b]
% to solve Rm(S,[a,b]) we need to minimize over t and maximize over u
uMode = 'max'; data0  = -data_original; a = 0.5, b =1;
tau_F = linspace(0,a,Nt);
schemeData.uMode = uMode;
[data_Rmfull, ~, ~] = ...
  HJIPDE_solve(data0, tau_F, schemeData, 'none', HJIextraArgs);
data_Rm = -squeeze(data_Rmfull(:,:,:,end));
figure(3)
    clf
    h0 = visSetIm(grid, data_Rmfull(:,:,:,1));
    h0.FaceAlpha = .3;    hold on
    h = visSetIm(grid, data_Rmfull(:,:,:,end));
    h.FaceAlpha = .3;
    title('The reachable set at the end')
    x0 = [0;0;0];
    if eval_u(grid,data_Rm,x0)>0
        plot3(0,0,0,'*r');
    else
        plot3(0,0,0,'*g');
    end
    hold off

% \bar{Rm(\bar{S},[a,b])}, which corresponds to G[a,b]S
% to solve Rm(S,[a,b]) we need to minimize over t and maximize over u

% method 1, first solve Rm(\bar{S},[a,b])
uMode = 'max'; data0 = -data_original; data_barS = -data0; a = 0.5, b =1;
tau_F = linspace(0,a,Nt);
schemeData.uMode = uMode;
[data_temp, ~, ~] = ...
  HJIPDE_solve(data_barS, tau_F, schemeData, 'none', HJIextraArgs);
data_RmBarS = -squeeze(data_temp(:,:,:,end)); % data_RmBarS>0 --> inside set

figure(4)
    clf
    h0 = visSetIm(grid, data_temp(:,:,:,1));
    h0.FaceAlpha = .3;    hold on
    h = visSetIm(grid, data_temp(:,:,:,end));
    h.FaceAlpha = .3;
    hold on
    title('The reachable set at the end')
    hold off
    
% then take the complement set 
data_G = -data_RmBarS;
figure(5)
    clf
    h0 = visSetIm(grid, data_original);
    h0.FaceAlpha = .3;    hold on
    h = visSetIm(grid, data_G);
    h.FaceAlpha = .3;
    hold on
    title('The set of G[a,b] S')
    x0 = [0;0;0];
    if eval_u(grid,data_G,x0)>0
        plot3(0,0,0,'*r');
    else
        plot3(0,0,0,'*g');
    end
    hold off
    
 % method 2, For unicycle model, G[a,b]S is also equivalent to F[0,a]S
uMode = 'min'; data0  = -data_original;
tau_F = linspace(0,a,Nt);
schemeData.uMode = uMode;
[data_F, ~, ~] = ...
  HJIPDE_solve(data0, tau_F, schemeData, 'zero', HJIextraArgs);
data_RM = -squeeze(data_F(:,:,:,end)); %tau

figure(6)
    clf
    h0 = visSetIm(grid, data_F(:,:,:,1));
    h0.FaceAlpha = .3;    hold on
    h = visSetIm(grid, data_F(:,:,:,end));
    h.FaceAlpha = .3;
    title('The reachable set at the end');
    x0 = [0;0;0];
    if eval_u(grid,data_RM,x0)>0
        plot3(0,0,0,'*r');
    else
        plot3(0,0,0,'*g');
    end
    hold off
    
    
 %% testing plot
 figure(7)
clf
h0 = visSetIm(grid, data0);
h0.FaceAlpha = .3;    hold on

h0 = visSetIm(grid, data1);
h0.FaceAlpha = .3; 
h0.FaceColor = 'cyan';

h0 = visSetIm(grid, data2);
h0.FaceAlpha = .3; 
h0.FaceColor = 'blue';
% h = visSetIm(grid, data_F(:,:,:,end));
% h.FaceAlpha = .3;
% title('The reachable set at the end');
% x0 = [0;0;0];
% if eval_u(grid,data_RM,x0)>0
%     plot3(0,0,0,'*r');
% else
%     plot3(0,0,0,'*g');
% end
hold off