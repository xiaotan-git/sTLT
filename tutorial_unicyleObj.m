clear all; close all;

%% build a unicycle CBF object

% Grid
grid_min = [-5; -5; -pi]; % Lower corner of computation domain
grid_max = [5; 5; pi];    % Upper corner of computation domain
N = [41; 41; 41];         % Number of grid points per dimension
pdDims = 3;               % 3rd dimension is periodic
g = createGrid(grid_min, grid_max, N, pdDims);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions

% Time interval
timeInterval=[1,2,3]; % timeInterval = [tb_low, tset_low, tb_upp]

% target set
R = 1; % radius
% data0 = shapeCylinder(grid,ignoreDims,center,radius)
% CBF construction for the set X = {x: data0(x)>=0} 
data0 = -shapeCylinder(g, 3, [0; 0; 0], R);
% also try shapeRectangleByCorners, shapeSphere, etc.

% input bounds
vRange = [-1 1];
wMax = 1;

% CBF construction for the set X = {x: data0(x)>=0} 
uniCBF = unicycleCBF(timeInterval,g,data0);

% test out the cbf value and its gradient
x = [1; 1; 0]; t = 1.2;
cbfValue = uniCBF.value(x,t)
cbfGrad = uniCBF.grad(x,t)


x = [0.4; 0.2; 0]; t = 2.2;
cbfValue = uniCBF.value(x,t)
cbfGrad = uniCBF.grad(x,t)

%% Do simulation

simulateTraj = true;
x = [0.8; 0.8; 0]; t = 1.2;
xint = x; tint = t;

enterflag = 0;

if simulateTraj
    
    initCBFValue = uniCBF.value(xint,tint);
%     if initCBFValue<0
%         error('initial CBF value is smaller than 0')
%     end


    dt = 0.01;
    MaxIter = floor((2.4 - tint)/dt);
    x_vec = zeros(3,MaxIter); x_vec(:,1) = xint;
    u_vec = zeros(2,MaxIter);
    t_vec = zeros(1,MaxIter); t_vec(1,1) = tint;
    b_vec = zeros(1,MaxIter);

    for iter = 1:MaxIter-1
        x = x_vec(:,iter);
        t = t_vec(1,iter);
        
        CBFval = uniCBF.value(x,t);
        CBFGrad = uniCBF.grad(x,t);
        
        fx = zeros(3,1); gx = [cos(x(3)) 0; sin(x(3)) 0; 0 1];
        A = CBFGrad(1:3,1)'*gx;
        b = 2*CBFval+CBFGrad(1:3,1)'*fx+CBFGrad(4,1);
        
        lb = -1.1*ones(2,1);
        ub = 1.1*ones(2,1);
        
        u = cbf_QP(A,b,lb,ub);        
        x_new = x+ dt*(fx+gx*u);
        t_new = t + dt;
        
        if x_new(3)>pi-1e-2
            x_new(3) = x_new(3) - 2*pi;
        else
            if x_new(3)<-pi+1e-2
                x_new(3) = 2*pi+x_new(3);
            end
        end
    
        x_vec(:,iter+1) = x_new;
        t_vec(:,iter+1) = t_new;
        u_vec(:,iter) = u;
        b_vec(:,iter) = CBFval;
        
        if norm(x_new(1:2,1))<1 && enterflag == 0
            disp(t)
            enterflag = 1;
        end
        
    end

    close all;
    figure(7)
    plot(x_vec(1,:),x_vec(2,:))
    title('trajectory and target set')
    hold on;
    [g2D, data2D] = proj(g, data0, [0 0 1]);
    visSetIm(g2D, data2D, 'green');

    figure(8)
    plot(t_vec,b_vec)
    title('barrier function value')
    
    figure(9)
    plot(t_vec,u_vec)
    title('control inputs')
    
end