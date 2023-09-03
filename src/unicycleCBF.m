classdef unicycleCBF
    %UNICYCLECBF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % current state in R^3
        x
        % current time in R+
        t
        
        % input bounds
        vRange
        wMax
        
        % time interval [tb_low, tset_low, tb_upp]
        timeInterval
        % grid over states
        grid
        g_with_time
        % target set 
        data0
        
        % value function (also CBF function)
        dataCBF
        
        % derivative function
        dataDeriv
        dataDerivPhase2

        % for sTLT tree spec
        temporalOperator
        setName
    end
    
    methods
        function obj = unicycleCBF(timeInterval,grid,data0,vRange,wMax)
            %UNICYCLECBF Construct an instance of this class
            if nargin == 3
                vRange = [-1 1];
                wMax = 1;
            end
            
            % To follow the convention in the Reachability toolbox
            % the target set X = {x: data0(x)<=0} 
            data0 = -data0; % careful, now we reverse the sign            
            % time grid
            tb_low = timeInterval(1); 
            tset_low = timeInterval(2);
            tb_upp =  timeInterval(3);
            
            %% time vector
            % if there exist two phrases
            if abs(tset_low - tb_low) >0.1
                Nt = 100;
                grid_min = grid.min; grid_max = grid.max; N = grid.N;
                pdDims = 3;
                g_with_time = createGrid([grid_min; tb_low], [grid_max;tset_low], [N; Nt], pdDims);

                t0 = 0;
                tMax = tset_low - tb_low; % todo: what if tset_low = tb_low
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
                HJIextraArgs.visualize.figNum = 100; %set figure number
                HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
                HJIextraArgs.quiet = true;

                % uncomment if you want to see a 2D slice
                %HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
                %HJIextraArgs.visualize.plotData.projpt = [0]; %project at theta = 0
                %HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

                % before calculating the CBF, let us do a re-initialization
                % more details, refer to Sec 2.2.1 The reinitialization equation in toolboxLS-1.1 manual

                %[data, tau, extraOuts] = ...
                % HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
                [data, ~, ~] = ...
                  HJIPDE_solve(data0, tau, schemeData, 'minVOverTime', HJIextraArgs);

                % now we flip the time and change the sign back
                dataCBF = -flip(data,4);
                dataDeriv = computeGradients(g_with_time, dataCBF);
                dataDerivPhase2 = computeGradients(grid, -data0); % reverse the sign here
                obj.g_with_time = g_with_time;
                obj.dataDerivPhase2 = dataDerivPhase2; % unique for 2 phases
            else
                % if only one phase exists
                dataCBF = - data0;
                dataDeriv = computeGradients(grid, dataCBF);
            end
            
            %% assign the computed parameters into the object
            obj.data0 = -data0; % reverse it back
            obj.dataCBF = dataCBF;
            obj.dataDeriv = dataDeriv;
            obj.grid = grid;
            obj.timeInterval = timeInterval;
            obj.vRange = vRange;
            obj.wMax = wMax;
            
        end
        
        function CBFval = value(obj,x,t)
            if nargin==1
                x = obj.x;
                t = obj.t;
            end
            
            if isequal(size(x),[3 1])
                x = x';
            end
            
            % if t is outside of the time domain
            if t<obj.timeInterval(1) || t>obj.timeInterval(3)
                error('time instant not within barrier timeInterval')
            end
            
            % if t is in phase 1 out of two phases
            if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1 && t>=obj.timeInterval(1) ...
                    && t<=obj.timeInterval(2)
                CBFval = eval_u(obj.g_with_time,obj.dataCBF,[x, t]);
            else
                % here we combine two cases (phase 2 or only 1 phase exists)
                CBFval = eval_u(obj.grid,obj.data0,x);
            end
                    
        end
        
        function CBFGrad = grad(obj,x,t)
            if nargin==1
                x = obj.x;
                t = obj.t;
            end
            
            if isequal(size(x),[3 1])
                x = x';
            end
            
            % if t is in phase 1 out of two phases
            if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1 && t>=obj.timeInterval(1) ...
                    && t<=obj.timeInterval(2)
                CBFGrad = eval_u(obj.g_with_time,obj.dataDeriv,[x, t]);
            else
                % if t is in phase 2 out of two phases
                if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1
                    CBFGradX = eval_u(obj.grid,obj.dataDerivPhase2,x);
                    CBFGrad = [CBFGradX; 0];
                else
                % if only one phase exists
                    CBFGradX = eval_u(obj.grid,obj.dataDeriv,x);
                    CBFGrad = [CBFGradX; 0];
                end
            end
            
        end

        function obj = updateDomain(obj,t,t_shift)
            % at time t, with time shift t_shift;
            if t<obj.timeInterval(3)
                new_timeInterval = obj.timeInterval - t_shift;
                % new_timeInterval(1) = obj.timeInterval(1);
                obj.timeInterval = new_timeInterval;
                
                % also update the time in g_with_time
                % check ToolboxLS/Kernel/Grids/processGrid.m for more details
                if ~isempty(obj.g_with_time)
                    gridIn.dim = obj.g_with_time.dim;
                    gridIn.min = obj.g_with_time.min - [0;0;0;t_shift];
                    gridIn.max = obj.g_with_time.max - [0;0;0;t_shift];
                    gridIn.N = obj.g_with_time.N;
                    gridIn.bdry = obj.g_with_time.bdry;
                    new_grid_with_time = processGrid(gridIn);

                    obj.g_with_time = new_grid_with_time;
                end

            end
            
        end

        function tf = isInTimeDomain(obj,t)
            if t>=obj.timeInterval(1) && t<=obj.timeInterval(3)
                tf = true;
            else
                tf = false;
            end
        end
    end
end

