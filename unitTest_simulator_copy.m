clear all; close all;
case_numb = 1;
switch case_numb
    case 1
        load('tree_singleIntegrator_example1.mat');
        load('cbf_cell_singleIntegrator_example1.mat');
    case 2
        load('tree_singleIntegrator_example2.mat');
        load('cbf_cell_singleIntegrator_example2.mat');
end



%alpha = [10 1 1 0.8 10 10 1 0.6];
% [1 1 1 1 10 10 1 1] 
% alpha = [10 1 1 1 10 10 1 1]; 
% alpha = [10 1 10 1 10 10 10 1]; % works for case 1
alpha = [1 10 1 1 1 1 10 1]; % works for case 2_1

%% define dynamics
dyn.f = @(x,t) [0;0]; 
dyn.g= @(x,t) [1 0; 0 1];


%% simulator
dt = 0.01;

if case_numb == 1
    T = 25;
end
if case_numb == 2
    T = 20;
end
t_vec = 0:dt:T;
total_iter = length(t_vec);
x_vec = zeros(2,total_iter);
% x_vec(:,1) = 10*rand(2,1)-5;
% x_vec(:,1) = [-6; 2]; 
x_vec(:,1) = [-2; 3.5]; % case 1
% x_vec(:,1) = [-1.7; 0.5]; % case 2
u_vec = zeros(2,total_iter);
unom_vec = zeros(2,total_iter);
I = diag(ones(2,1));

A_vec = zeros(length(cbf_cell),2,total_iter);
b_vec_log = zeros(length(cbf_cell),total_iter);

bFunc_vec = zeros(length(cbf_cell),total_iter);

for iter = 1:total_iter-1
    x = x_vec(:,iter);
    t = t_vec(:,iter);
    % stack all the CBF conditions. Here we follow Au+b >=0
    A = zeros(length(cbf_cell),2); b_vec = zeros(length(cbf_cell),1);
    bFunc =  zeros(length(cbf_cell),1);
    fx = dyn.f(x,t); gx = dyn.g(x,t);
    
    % check for possible updates of ts
    for i = 1:length(tree.temporalFragments)
        setNodeIndex2Check = tree.temporalFragments{i}(2);
        setNode2Check = tree.nodeList{setNodeIndex2Check};
        if setNode2Check.isInStartInterval(t) && setNode2Check.isInRegion(x) ...
                && (setNode2Check.startTimeInterval(2)-setNode2Check.startTimeInterval(1)>=1e-3)
            disp('Now at time '); disp(t); disp(' update the set node')
            disp(setNode2Check.nodeName)
            % record the bar{t}_s(Xj)
            t1 = setNode2Check.startTimeInterval(2);
            
            setNodeIndex2Update = setNodeIndex2Check;
            % Algorithm 4: online update the start time interval 
            set2UpdateIndexList = tree.updateStartTimeInterval(setNodeIndex2Update,t);
            
            % update CBFs online
            disp('CBFs to be updated are')
            for i_set = 1:length(set2UpdateIndexList)
                for i_cbf = 1:length(cbf_cell)
                    % check cbf index to update 
                    setNameCBF = cbf_cell{i_cbf}.setName;
                    setIndexFromCBF = tree.findSetNode(setNameCBF);
                    if set2UpdateIndexList{i_set} == setIndexFromCBF 
                        % here we use setNodeName as an identifier for temporal fragments
                        disp(i_cbf)
                        cbf_cell{i_cbf} = cbf_cell{i_cbf}.updateDomain(t,t1-t); % at time t, with shift t1 - t;
                    end
                end
            end
        end
    end
    
    
    
    
    for i = 1:length(cbf_cell)
        tmin = cbf_cell{i}.timeInterval(1);
        tmax = cbf_cell{i}.timeInterval(3);
        if t >= tmin && t<= tmax
            bFunc(i) = cbf_cell{i}.value(x,t);
            bgrad = cbf_cell{i}.grad(x,t);
            A(i,:) = bgrad(1:2,1)'*gx;
            b_vec(i)= alpha(i)*bFunc(i) + bgrad(1:2,1)'*fx+bgrad(3,1);
            
        end
    end
    
     lb = [];up = [];
%     lb = -1*ones(2,1); up = 1*ones(2,1);
     
    % heuristical nominal control 
    % Intuitive 1: try generate unom that fulfill all cbfs
%     A_copy = A;
%     ind = (vecnorm(A_copy')<=1e-2)';
%     A_copy(ind,:)=[];
%     if rank(A_copy)== size(A_copy,1)
%         unom = 2*A_copy'*inv(A_copy*A_copy')*ones(size(A_copy,1),1);
%     else
%         unom = zeros(size(A_copy,2),1);
%     end

    % Intuitive 2: try generate unom when only one cbf exists
%     if rank(A_copy)== 1
%         unom = A_copy';
%     else
%         unom = zeros(size(A_copy,2),1);
%     end

    % Intuitive 3: generate unom when only one F(a,b)X is active
%     activeCBFList = {};
%     for i = 1:length(cbf_cell)
%         if cbf_cell{i}.isInTimeDomain(t)
%             activeCBFList{end+1} = cbf_cell{i};
%         end
%     end
%     %check how many F(a,b)X is active
%     F_num = 0;activeFIndex = 0;
%     if numel(activeCBFList) >0
%         for i = 1:length(activeCBFList)
%             temporalOperator = activeCBFList{i}.temporalOperator;
%             if strcmp(temporalOperator(1),'F')
%                 F_num = F_num +1;
%                 activeFIndex = i;
%             end
%         end
%     end
%     if F_num==1
%         activeFCBF = activeCBFList{i};
%         bgrad = activeFCBF.grad(x,t);
%         ai = bgrad(1:2,1)'*gx;
%         unom = 0.1*ai';
%     else
%         unom = zeros(size(A_copy,2),1);
%     end
    
    % intuitive 4: guide the trajectory towards fulfilling the CBF with the nearest ddl
    % if two sets with the same ddl, prioritize the CBF corresponding to an F operator
    nearestDDL = inf; nearestDDLInd = 0; nearestDDLType = '';
    for i = 1:length(cbf_cell)
        if cbf_cell{i}.isInTimeDomain(t) && (nearestDDL>=cbf_cell{i}.timeInterval(end))
            %  CBF corresponding to an F operator has priority
            if cbf_cell{i}.temporalOperator(1) == 'F' || (cbf_cell{i}.temporalOperator(1) == 'G' && ...
                isempty(nearestDDLType))
                nearestDDL = cbf_cell{i}.timeInterval(end);
                nearestDDLInd = i;
                nearestDDLType = cbf_cell{i}.temporalOperator(1);
            end
        end
    end
    if nearestDDLInd ~= 0
        nearesrDDLCBF = cbf_cell{nearestDDLInd};
        bgrad = nearesrDDLCBF.grad(x,t);
        ai = bgrad(1:2,1)'*gx;
        unom = 0.1*ai';
    else
        unom = zeros(size(A,2),1);
    end
    
%     unom = zeros(size(A,2),1);
    
    % cal u
    try
        u = cbf_QP(A,b_vec,lb,up,unom);
    catch
        u = unom;
        disp('QP not feasible at time')
        disp(iter*dt);
    end
        
    % hard saturation
    lb = -1*ones(2,1); up = 1*ones(2,1);
    if ~(all(u>=lb) && all(u<=up))
        u = sat_func(u,lb,up); 
%         disp('u from QP not feasible. Saturation operation is done at time')
%         disp(iter*dt);
    end
    % logging
    x_vec(:,iter+1) = x_vec(:,iter)+dt*(fx+gx*u);
    u_vec(:,iter) = u;
    unom_vec(:,iter) = unom;
    A_vec(:,:,iter) = A;
    b_vec_log(:,iter) = b_vec;
    bFunc_vec(:,iter) = bFunc;
end

%% plotting
figure
hold on;
for i = 1:numel(cbf_cell)
    cbf_i = cbf_cell{i};
    drawCircle(cbf_i.circularCenter, cbf_i.circularRadius);
    st= [cbf_i.setName(1:2) '{' cbf_i.setName(3:end) '}'];
    text(cbf_i.circularCenter(1,1)+ cbf_i.circularRadius -1,cbf_i.circularCenter(2,1),st,...
    'VerticalAlignment','middle',...                
    'HorizontalAlignment', 'left',... % alignment with the node
    'FontSize', 12);
end
plot(x_vec(1,1),x_vec(2,1),'ko','MarkerSize',12);
% plot3(x_vec(1,:),x_vec(2,:),t_vec,'linewidth',2);

% to plot with a time-varying color
% for i = 1: size(x_vec_cell,2)
%     x_vec = x_vec_cell{i};
%     plot(x_vec(1,1),x_vec(2,1),'ko','MarkerSize',12);
%     plot_traj_color(x_vec,t_vec);
% end
% colorbar('Ticks',[0,5,10,15,20],...
%          'TickLabels',{'0s','5s','10s','15s','20s'})
% colormap('turbo')
plot_traj_color(x_vec,t_vec);
colorbar('Ticks',[0,5,10,15,20],...
         'TickLabels',{'0s','5s','10s','15s','20s'})
colormap('turbo')
sample = find(t_vec==1)-1;
t_text = t_vec(1:sample:end);
st_text_pre = num2str(t_text'); st_text_app = st_text_pre(:,1); st_text_app(:,1) = 's';
st_text = [st_text_pre st_text_app];
plot(x_vec(1,1:sample:end),x_vec(2,1:sample:end),'square')
text(x_vec(1,1:sample:end)+0.05,x_vec(2,1:sample:end),st_text,...
    'VerticalAlignment','middle',...                
    'HorizontalAlignment', 'left',... % alignment with the node
    'FontSize', 8);
axis equal


figure
plot(t_vec,x_vec(1,:));
hold on;
plot(t_vec,x_vec(2,:));
legend('x1','x2');
title("x(t)")

figure
plot(t_vec,u_vec)
title("u(t)")


figure
plot(t_vec,unom_vec)
title("unom(t)")

figure
plot(t_vec,bFunc_vec)
legend;
title("barriers(t)")