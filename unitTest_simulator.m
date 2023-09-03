% clear all; close all;
case_numb = 2;
dyn_model_name = 'singleIntegrator'; % 'singleIntegrator' or 'unicycle'
switch case_numb
    case 1
        load(['tree_' dyn_model_name '_example1.mat']);
        load(['cbf_cell_' dyn_model_name '_example1.mat']);
    case 2
        load(['tree_' dyn_model_name '_example2.mat']);
        load(['cbf_cell_' dyn_model_name '_example2.mat']);
end
disp(['Now running the case ' num2str(case_numb) ...
    ' of the model ' dyn_model_name])


%alpha = [10 1 1 0.8 10 10 1 0.6];
% [1 1 1 1 10 10 1 1] 
% alpha = [10 1 1 1 10 10 1 1]; 
% alpha = [10 1 10 1 10 10 10 1]; % works for case 1
alpha = [1 10 1 1 1 1 10 1 1]; % works for case 2_1

%% dynamics related parameters
if strcmp(dyn_model_name,'singleIntegrator')
    dyn.f = @(x,t) [0;0]; 
    dyn.g= @(x,t) [1 0; 0 1];
    size_state = 2;
    size_input = 2;
    uMax = 1;
    % lb = -uMax*ones(size_input,1); up = uMax*ones(size_input,1);

    % case 1
%     x0 = [-2; 3.5]; % [-6; 2] or [-2; 3.5] or [-1.7; 0.5] or [-20; -5];
    
    % case 2
    x0 = [-1.7; 0.5]; % [-2; 3.5] or [-1.7; 0.5]
    
    
    gain_nom = 0.1;
end
if strcmp(dyn_model_name,'unicycle')
    dyn.f = @(x,t) [0;0;0]; 
    dyn.g= @(x,t) [cos(x(3)) 0; sin(x(3)) 0; 0 1];
    size_state = 3;
    size_input = 2;
    wMax = 1; vMax = 1;
    lb = -[vMax;wMax]; up = [vMax;wMax];

    x0 =[-2; 3.5;pi/2]; % [-2; 3.5;pi/2] or [-6;2;0]
    gain_nom = 0.01; % a tuning parameter
end


%% simulator
dt = 0.01;

if case_numb == 0
    T = 15;
end

if case_numb == 1
    T = 25;
end
if case_numb == 2
    T = 20;
end
t_vec = 0:dt:T;
total_iter = length(t_vec);
x_vec = zeros(size_state,total_iter);
% x_vec(:,1) = 10*rand(2,1)-5;
% x_vec(:,1) = [-6; 2]; 
x_vec(:,1) = x0; % case 1
% x_vec(:,1) = [-1.7; 0.5]; % case 2
u_vec = zeros(size_input,total_iter);
unom_vec = zeros(size_input,total_iter);

A_vec = zeros(length(cbf_cell),size_input,total_iter);
b_vec_log = zeros(length(cbf_cell),total_iter);

bFunc_vec = zeros(length(cbf_cell),total_iter);

% for iter = 1:2
for iter = 1:total_iter-1
    x = x_vec(:,iter);
    t = t_vec(:,iter);
    % stack all the CBF conditions. Here we follow Au+b >=0
    A = zeros(length(cbf_cell),2); b_vec = zeros(length(cbf_cell),1);
    bFunc =  zeros(length(cbf_cell),1);
    fx = dyn.f(x,t); gx = dyn.g(x,t);
    
    % check for possible updates of ts
    for i = 1:length(cbf_cell)
        setName = cbf_cell{i}.setName;
        setNodeIndex2Check = tree.findSetNode(setName);
        setNode2Check = tree.nodeList{setNodeIndex2Check};
        if setNode2Check.isInStartInterval(t) && setNode2Check.isInRegion(x) ...
                && (setNode2Check.startTimeInterval(2)-setNode2Check.startTimeInterval(1)>=1e-3)
            disp(['Now at time ' num2str(t)  's, update the set node ' setNode2Check.nodeName])

            % record the bar{t}_s(Xj)
            t1 = setNode2Check.startTimeInterval(2);
            
            setNodeIndex2Update = setNodeIndex2Check;
            % Algorithm 4: online update the start time interval 
            set2UpdateIndexList = tree.updateStartTimeInterval(setNodeIndex2Update,t);
            
            % update CBFs online
            for i_set = 1:length(set2UpdateIndexList)
                for i_cbf = 1:length(cbf_cell)
                    % check cbf index to update 
                    setNameCBF = cbf_cell{i_cbf}.setName;
                    setIndexFromCBF = tree.findSetNode(setNameCBF);
                    if set2UpdateIndexList{i_set} == setIndexFromCBF 
                        % here we use setNodeName as an identifier for temporal fragments
                        cbf_cell{i_cbf} = cbf_cell{i_cbf}.updateDomain(t,t1-t); % at time t, with shift t1 - t;
                        disp(['CBFs to be updated are ' ...
                            cbf_cell{i_cbf}.temporalOperator cbf_cell{i_cbf}.setName ])
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
            A(i,:) = bgrad(1:end-1,1)'*gx;
            b_vec(i)= alpha(i)*bFunc(i) + bgrad(1:end-1,1)'*fx+bgrad(end,1);
            
        end
    end
    
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
    
    % intuition 4: guide the trajectory towards fulfilling the CBF with the nearest ddl
    % if two sets with the same ddl, prioritize the CBF corresponding to an F operator
    nearestDDL = inf; nearestDDLInd = 0; nearestDDLType = '';
    % for i = 1:length(cbf_cell)
    %     if cbf_cell{i}.isInTimeDomain(t) && (nearestDDL>cbf_cell{i}.timeInterval(end))
    %         % if cbf_i has the nearest DDL
    %         nearestDDL = cbf_cell{i}.timeInterval(end);
    %         nearestDDLInd = i;
    %         nearestDDLType = cbf_cell{i}.temporalOperator(1);
    %     else
    %         % if there is a tie in DDL
    %         if cbf_cell{i}.isInTimeDomain(t) && (nearestDDL==cbf_cell{i}.timeInterval(end))
    %             %  CBF corresponding to an F operator has priority
    %             if cbf_cell{i}.temporalOperator(1) == 'F' || (cbf_cell{i}.temporalOperator(1) == 'G' && ...
    %                 isempty(nearestDDLType))
    %                 nearestDDL = cbf_cell{i}.timeInterval(end);
    %                 nearestDDLInd = i;
    %                 nearestDDLType = cbf_cell{i}.temporalOperator(1);
    %             end
    %         end
    %     end
    % end
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
        ai = bgrad(1:end-1,1)'*gx;
        unom = ai'/(norm(ai)+gain_nom); 
    else
        unom = zeros(size_input,1);
    end
    
    %   unom = zeros(size(A,2),1);
    
    % cal u
    lb = []; up = [];
    try
        u = cbf_QP(A,b_vec,lb,up,unom);
    catch
        u = unom;
        disp(['QP not feasible at time ' num2str(iter*dt)])
        disp(['Current iteration: ' num2str(iter)])
    end
        
    % hard saturation
    lb = -1.0*ones(size_input,1); up = 1.0*ones(size_input,1);
    if ~(all(u>=lb) && all(u<=up))
        u = sat_func(u,lb,up); 
        % disp('u from QP not feasible. Saturation operation is done at time')
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
leafNodeIndex = tree.findLeafNode();
for i = 1:numel(cbf_cell)
    cbf_i = cbf_cell{i};
    % only plot the predicate region
    setNodeIndex = tree.findSetNode(cbf_i.setName);
    if any(leafNodeIndex==setNodeIndex)
        if any(strcmp(properties(cbf_i), 'circularCenter'))
            drawCircle(cbf_i.circularCenter, cbf_i.circularRadius);
            st= [cbf_i.setName(1:2) '{' cbf_i.setName(3:end) '}'];
            text(cbf_i.circularCenter(1,1)+ cbf_i.circularRadius -1,cbf_i.circularCenter(2,1),st,...
                'VerticalAlignment','middle',...                
                'HorizontalAlignment', 'left',... % alignment with the node
                'FontSize', 12);
        end

        if any(strcmp(properties(cbf_i), 'data0'))
            [g2D, data2D] = proj(cbf_i.grid, cbf_i.data0, [0 0 1]);
            visSetIm(g2D, data2D, 'green');
        end
    end
end
plot(x_vec(1,1),x_vec(2,1),'ko','MarkerSize',12);
% plot3(x_vec(1,:),x_vec(2,:),t_vec,'linewidth',2);



plot(x_vec(1,:),x_vec(2,:));
sample = find(t_vec==1)-1;
t_text = t_vec(1:sample:end);
st_text_pre = num2str(t_text'); st_text_app = st_text_pre(:,1); st_text_app(:,1) = 's';
st_text = [st_text_pre st_text_app];
plot(x_vec(1,1:sample:end),x_vec(2,1:sample:end),'x')
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

%% plot the figures in the paper
% data logging
% x_vec_cell = {};
% x_vec_cell{end+1} = x_vec;

% to plot with a time-varying color
figure
for i = 1: size(x_vec_cell,2)
    x_vec = x_vec_cell{i};
    plot(x_vec(1,1),x_vec(2,1),'ko','MarkerSize',12);
    plot_traj_color(x_vec,t_vec);
end
colorbar('Ticks',[0,5,10,15,20],...
         'TickLabels',{'0s','5s','10s','15s','20s'})
colormap('turbo')



 % adjusting the figure
if case_numb == 1
    axis equal
    xlim([-8,8])
    ylim([-7,4])
    box on;
    ax = gca;
    ax.FontSize = 12; 
    xlabel('$x_1$','interpreter','latex','fontsize',16)
    ylabel('$x_2$','interpreter','latex','fontsize',16)
    
    X5.c = [-4;-4]; X5.r = 1; 
    X3.c = X5.c; X3.r = 3;
    % X1.c = X5.c; X1.r = 18;
    % drawCircle(X1.c, X1.r);
    X9.c = [1;-4]; X9.r = 2;
    X8.c = [4;0]; X8.r = 4;
    % drawCircle(X8.c, X8.r);
    % drawCircle(X9.c, X9.r);
    X6.c = [4;0]; X6.r = 4;
    X7.c = [1;-4]; X7.r = 12;
    X4.c = [4;0]; X4.r = 4;
    X2.c = [4;0]; X2.r = 19;
    
    % drawCircle(X3.c, X3.r);
    drawCircle(X4.c, X4.r);
    drawCircle(X5.c, X5.r);
    drawCircle(X8.c, X8.r);
    drawCircle(X9.c, X9.r);
    drawCircle(X3.c, X3.r);
    
    delete(findall(gcf,'type','annotation'))
    set(gcf,'Units','normalized');
    [xf5,yf5] = axescoord2figurecoord([X5.c(1)-0.2;X5.c(1)+0.2],[X5.c(2)-0.2;X5.c(2)+0.2]);
    annotation('textbox',[xf5(1)-0.02 yf5(1) xf5(2)-xf5(1) yf5(2)-yf5(1)],'interpreter','latex','String','$X_{5}$','EdgeColor','None','fontsize',16)
    [xf4,yf4] = axescoord2figurecoord([X4.c(1)-0.2;X4.c(1)+0.2],[X4.c(2)-0.2;X4.c(2)+0.2]);
    annotation('textbox',[xf4(1) yf4(1) xf4(2)-xf4(1) yf4(2)-yf4(1)],'interpreter','latex','String','$X_{8} / X_{4}$','EdgeColor','None','fontsize',16)
    [xf3,yf3] = axescoord2figurecoord([X3.c(1)-0.15;X3.c(1)+0.2],[X3.c(2)-0.2;X3.c(2)+0.2]);
    annotation('textbox',[xf3(1)-0.11 yf3(1) xf3(2)-xf3(1) yf3(2)-yf3(1)],'interpreter','latex','String','$X_{3}$','EdgeColor','None','fontsize',16)
    [xf9,yf9] = axescoord2figurecoord([X9.c(1)-0.2;X9.c(1)+0.2],[X9.c(2)-0.2;X9.c(2)+0.2]);
    annotation('textbox',[xf9(1) yf9(1) xf9(2)-xf9(1) yf9(2)-yf9(1)],'interpreter','latex','String','$X_{9}$','EdgeColor','None','fontsize',16)
end

if case_numb == 2
    axis equal
    xlim([-4,4])
    ylim([-2,4])
    box on;
    ax = gca;
    ax.FontSize = 12; 
    xlabel('$x_1$','interpreter','latex','fontsize',16)
    ylabel('$x_2$','interpreter','latex','fontsize',16)
    
    X1.c = [-2;2]; X1.r = 1; 
    % X2.c = [2;2]; X2.r = 1;
    X3.c = [-1;0];X3.r = 1.5;
    X4.c = [0;1.5]; X4.r = 1;
    X2.c = [1.5;0.5]; X2.r = 1;
    Square5.c = [-3; -1.5]; Square5.r = 0.5;
    Square5.obstacle = false;Square5.rep = 'square';
    Obs = X4;
    
    drawCircle(X1.c, X1.r);
    drawCircle(X2.c, X2.r);
    % drawCircle(X22.c, X22.r);
    drawCircle(X3.c, X3.r);
    drawCircle(X4.c, X4.r,1);
    drawBox(Square5.c,Square5.r);

    delete(findall(gcf,'type','annotation'))
    set(gcf,'Units','normalized');
    [xf1,yf1] = axescoord2figurecoord([X1.c(1)-0.2;X1.c(1)+0.2],[X1.c(2)-0.2;X1.c(2)+0.2]);
    annotation('textbox',[xf1(1)-0.02 yf1(1) xf1(2)-xf1(1) yf1(2)-yf1(1)],'interpreter','latex','String','$S_{\mu_1}$','EdgeColor','None','fontsize',16)
    [xf2,yf2] = axescoord2figurecoord([X2.c(1)-0.2;X2.c(1)+0.2],[X2.c(2)-0.2;X2.c(2)+0.2]);
    annotation('textbox',[xf2(1) yf2(1) xf2(2)-xf2(1) yf2(2)-yf2(1)],'interpreter','latex','String','$S_{\mu_2}$','EdgeColor','None','fontsize',16)
    [xf3,yf3] = axescoord2figurecoord([X3.c(1)-0.15;X3.c(1)+0.2],[X3.c(2)-0.2;X3.c(2)+0.2]);
    annotation('textbox',[xf3(1)-0.11 yf3(1) xf3(2)-xf3(1) yf3(2)-yf3(1)],'interpreter','latex','String','$S_{\mu_3}$','EdgeColor','None','fontsize',16)
    [xf4,yf4] = axescoord2figurecoord([X4.c(1)-0.2;X4.c(1)+0.2],[X4.c(2)-0.2;X4.c(2)+0.2]);
    annotation('textbox',[xf4(1) yf4(1) xf4(2)-xf4(1) yf4(2)-yf4(1)],'interpreter','latex','String','$S_{\mu_4}$','EdgeColor','None','fontsize',16)
    [xf5,yf5] = axescoord2figurecoord([Square5.c(1)-0.2;Square5.c(1)+0.2],[Square5.c(2)-0.2;Square5.c(2)+0.2]);
    annotation('textbox',[xf5(1) yf5(1) xf5(2)-xf5(1) yf5(2)-yf5(1)],'interpreter','latex','String','$S_{\mu_5}$','EdgeColor','None','fontsize',16)
end

% to plot with explicit time-stamps
timeStamp_interval = [2,20;
    4,20;
    6,18;
    8,18];

for i = 1:numel(x_vec_cell)
    x_vec = x_vec_cell{i};
    % stamp every 2s
    sample = find(t_vec==1)-1; 
%     start_iter = timeStamp_interval(i,1)*sample/2;
%     termin_iter = timeStamp_interval(i,2)*sample/2;
    t_text = [timeStamp_interval(i,1):2:timeStamp_interval(i,2)];
    st_text_pre = num2str(t_text'); st_text_app = st_text_pre(:,1); st_text_app(:,1) = 's';
    st_text = [st_text_pre st_text_app];
    for j = 1:numel(t_text)
        iter_timeStamp = sample*t_text(j);
        plot(x_vec(1,iter_timeStamp),x_vec(2,iter_timeStamp),'x')
        text(x_vec(1,iter_timeStamp)+0.05,x_vec(2,iter_timeStamp),st_text(j,:),...
            'VerticalAlignment','middle',...                
            'HorizontalAlignment', 'left',... % alignment with the node
            'FontSize', 8);
    end
end

for i = 1:2
    x_vec = x_vec_cell{i};
    text(x_vec(1,1)+0.05,x_vec(2,1),'0s',...
        'VerticalAlignment','middle',...                
        'HorizontalAlignment', 'left',... % alignment with the node
        'FontSize', 8);
end

% % if we want to have a subfigure
% % ax2 = axes('Position',[0.66 0.25 0.18 0.12]);
% ax2 = axes('Position',[0.08 0.75 0.18 0.06]);
% 
% x_vec = x_vec_cell{3};
% plot(x_vec(1,1),x_vec(2,1),'ko','MarkerSize',8);
% plot_traj_color(x_vec,t_vec);
% axis equal
% xlim([-6.1,-5.8]);
% ylim([1.85,2.05])
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% 
% 
% % timeStamp_interval = [4,7];
% % stamp every 2s
% sample = find(t_vec==1)-1; 
% %     start_iter = timeStamp_interval(i,1)*sample/2;
% %     termin_iter = timeStamp_interval(i,2)*sample/2;
% t_text = [timeStamp_interval(1,1):3:timeStamp_interval(1,2)];
% st_text_pre = num2str(t_text'); st_text_app = st_text_pre(:,1); st_text_app(:,1) = 's';
% st_text = [st_text_pre st_text_app];
% for j = 1:numel(t_text)
%     iter_timeStamp = sample*t_text(j);
%     plot(x_vec(1,iter_timeStamp),x_vec(2,iter_timeStamp),'x')
%     text(x_vec(1,iter_timeStamp)+0.05,x_vec(2,iter_timeStamp),st_text(j,:),...
%         'VerticalAlignment','middle',...                
%         'HorizontalAlignment', 'left',... % alignment with the node
%         'FontSize', 8);
% end
