function [grid,data] = uni_ReachableSetOperation(operation,time_interval,grid,data_original,dyn)
%UNI_REACHABLESETOPERATION This function takes in the set operation, and
%return the reachable set.
%   The original set is taken as the super-level set of data 0



    
% Put grid and dynamic systems into schemeData
schemeData.grid = grid;
schemeData.dynSys = dyn;
schemeData.accuracy = 'high'; %set accuracy


%% Compute value function

%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 1;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.quiet = true;


if operation == 'G'
    disp('calculating the set S1 by solving G[a,b] S0');
    [grid,data] = G_operation(time_interval,grid,data_original,schemeData,HJIextraArgs);
else
    if operation == 'F'
        disp('calculating the set S1 by solving F[a,b] S0');
        [grid,data] = F_operation(time_interval,grid,data_original,schemeData,HJIextraArgs);
    end
end

end

function [grid,data_RM] = G_operation(time_interval,grid,data_original,schemeData,HJIextraArgs)

% method 2, For unicycle model, G[a,b]S is also equivalent to F[0,a]S
uMode = 'min'; data0  = -data_original; Nt = 100;
a = time_interval(1);
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

end


function [grid,data_RM] = F_operation(time_interval,grid,data_original,schemeData,HJIextraArgs)
% RM(S,[a,b]), which corresponds to F[a,b]S
% to solve RM(S,[a,b]) we need to minimize over t and minimize over u
uMode = 'min'; data0  = -data_original; a = time_interval(1); b =time_interval(2); Nt = 100;
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

end


