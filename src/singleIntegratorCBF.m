classdef singleIntegratorCBF
    %SINGLEINTEGRATORCBF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % current state in R^2
        x
        % current time in R+
        t
        
        % input bounds
        vMax
        
        % time interval [tb_low, tset_low, tb_upp]
        timeInterval
        % target set when t = tset_low 
        circularCenter
        circularRadius
        
        ObsFlag
        
        % corresponding temporal fragment
        temporalOperator % 'F[a,b]'
        setName % 'X3'
    end
    
    methods
        function obj = singleIntegratorCBF(timeInterval,c,r,vMax,obs)
            %SINGLEINTEGRATORCBF Construct an instance of this class
            %   Detailed explanation goes here
            if nargin ==5
                if strcmp(obs,'obs') || (obs == true)
                    obj.ObsFlag = 1;
                else
                    obj.ObsFlag = 0;
                end
            else
                obj.ObsFlag = 0;
            end
            
            if nargin <4
                vMax = 1;
            end
            
            obj.vMax = vMax;
            obj.timeInterval = timeInterval;
            obj.circularCenter = c;
            obj.circularRadius = r;
        end
        
        function CBFval = value(obj,x,t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if nargin==1
                x = obj.x;
                t = obj.t;
            end
            
            if t<obj.timeInterval(1) || t>obj.timeInterval(3)
                error('time instant not within barrier timeInterval')
            end
            
            % if t is in phase 1 out of two phases
            if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1 ...
                    && t>=obj.timeInterval(1) && t<=obj.timeInterval(2)
                T = obj.timeInterval(2) - obj.timeInterval(1);
                CBFval = (obj.circularRadius+obj.vMax*T -obj.vMax*(t-obj.timeInterval(1)))^2 ...
                    - (x-obj.circularCenter)'*(x-obj.circularCenter);
            else
                % here we combine two cases (phase 2 or only 1 phase exists)
                CBFval = obj.circularRadius^2 ...
                    - (x-obj.circularCenter)'*(x-obj.circularCenter);
            end
            
            if obj.ObsFlag == 1
                CBFval =  - CBFval;
            end

        end
        
        function CBFGrad = grad(obj,x,t)
            if nargin==1
                x = obj.x;
                t = obj.t;
            end
            
            if t<obj.timeInterval(1) || t>obj.timeInterval(3)
                error('time instant not within barrier timeInterval')
            end
            
            % if t is in phase 1 out of two phases
            if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1 ...
                    && t>=obj.timeInterval(1) && t<=obj.timeInterval(2)
                T = obj.timeInterval(2) - obj.timeInterval(1);
                CBFGradX = -2*(x-obj.circularCenter);
                CBFGradt = 2*(obj.circularRadius+obj.vMax*T ...
                    -obj.vMax*(t-obj.timeInterval(1)))*(-obj.vMax);
                CBFGrad = [CBFGradX; CBFGradt];
            else
                % here we combine two cases (phase 2 or only 1 phase exists)
                CBFGradX = -2*(x-obj.circularCenter);
                CBFGradt = 0;
                CBFGrad = [CBFGradX; CBFGradt];
            end
            
            if obj.ObsFlag == 1
                CBFGrad =  - CBFGrad;
            end
            
        end
        
        function obj = updateDomain(obj,t,t_shift)
            % at time t, with time shift t_shift;
            if t<obj.timeInterval(3)
                new_timeInterval = obj.timeInterval - t_shift;
                % new_timeInterval(1) = obj.timeInterval(1);
                obj.timeInterval = new_timeInterval;
            end
            
        end
        
        function tf = isInTimeDomain(obj,t)
            if t>=obj.timeInterval(1) && t<=obj.timeInterval(3)
                tf = true;
            else
                tf = false;
            end
        end

        function region = regionAtTimeT(obj,t)
            % output the superlevel set region at time t
            if t<obj.timeInterval(1) || t>obj.timeInterval(3)
                error('time instant not within barrier timeInterval')
            end
            
            if obj.ObsFlag == 0
                % not an obstacle region
                % if t is in phase 1 out of two phases
                if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1 ...
                        && t>=obj.timeInterval(1) && t<=obj.timeInterval(2)
                    T = obj.timeInterval(2) - t;
                    region.c = obj.circularCenter;
                    region.r = obj.circularRadius+obj.vMax*T;
                else
                    % here we combine two cases (phase 2 or only 1 phase exists)
                    region.c = obj.circularCenter;
                    region.r = obj.circularRadius;
                end
                region.obstacle = obj.ObsFlag;
            else
                % if t is in phase 1 out of two phases
                if abs(obj.timeInterval(2) - obj.timeInterval(1)) >0.1 ...
                    && t>=obj.timeInterval(1) && t<=obj.timeInterval(2)
                    T = obj.timeInterval(2) - t;
                    region.c = obj.circularCenter;
                    region.r = max(obj.circularRadius-obj.vMax*T,0);
                else
                    % here we combine two cases (phase 2 or only 1 phase exists)
                    region.c = obj.circularCenter;
                    region.r = obj.circularRadius;
                end
                region.obstacle = obj.ObsFlag;
            end

            
        end
    end
end

