function new_region = regionCal(region,operatorSt,varargin)
%REGIONCAL This function calculates the region in the setNode that is the
%parent of the operatorNode (operator) and the setNode (region)
%legal operators inclde 'F[0,1]', 'G[0,1]', 'Vee', 'Wedge'
%for now only circular region for singleIntegrator dynamics is implemented
% umax = 1 by assumption
% possible input forms
%  new_region = regionCal(region1,operatorSt,region2,umax,regionRep)
%  e.g., new_region = regionCal(region,'F[0,1]')
%  new_region = regionCal(region1,'Vee',region2)
%  new_region = regionCal(region1,'Vee',region2,'umax',2)
%  new_region = regionCal(region1,'Vee',region2,'umax',2,'regionRep','circluar');


p = inputParser;
% % Define the expected parameters and their default values
% more details, see https://se.mathworks.com/help/matlab/ref/inputparser.html
% addParameter(p, 'paramName', defaultValue, validationFunction);
% optional propositional argument
addOptional(p, 'region2', '', @(x) isstruct(x));
% optional name-value pair argument
addParameter(p, 'umax', 1, @(x) isnumeric(x));
addParameter(p, 'vRange', [-1 1], @(x) isnumeric(x));
addParameter(p, 'wMax', 1, @(x) isnumeric(x));


% possible regionRep: circular, grid-data, levelset 
addParameter(p,'regionRep','circular',@(x) ischar(x));


parse(p, varargin{:});
% % parsedValue = p.Results.paramName;
% region = p.Results.region1;
% operatorSt = p.Results.operatorSt;
region2 = p.Results.region2;
umax = p.Results.umax;
regionRep = p.Results.regionRep;
vRange = p.Results.vRange;
wMax = p.Results.wMax;

switch regionRep
    case 'circular'
        new_region.c = [];
        new_region.r = [];
        new_region.obstacle = [];

        if operatorSt(1) == 'F'
            interval = str2num(operatorSt(2:end));
            new_region = region;

            if ~new_region.obstacle
                new_region.r = new_region.r + umax*interval(end);
            else
                new_region.r = max(0,new_region.r - umax*interval(end));
            end
                
        end

        if operatorSt(1) == 'G'
            interval = str2num(operatorSt(2:end));
            new_region = region;
            
            if ~region.obstacle
                new_region.r = new_region.r + umax*interval(1);
            else
                new_region.r = max(0,new_region.r - umax*interval(1));
            end
                
        end

        if operatorSt(1) == 'V'
            
            new_region = setUnion(region,region2);
                    
        end

        if operatorSt(1) == 'W'
            
            new_region = setIntersection(region,region2);
                    
        end
    case 'grid-data'
        % for backward reachability calculation
        t0 = 0; Nt = 41;

        %% reachability problem parameters
        uMode = 'min';
        dUni = Plane([0,0,0],wMax,vRange);

        % Put grid and dynamic systems into schemeData
        schemeData.grid = region.grid;
        schemeData.dynSys = dUni;
        schemeData.accuracy = 'high'; %set accuracy
        schemeData.uMode = uMode;

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

        if operatorSt(1) == 'F'
            new_region = region;
            interval = str2num(operatorSt(2:end));
            tMax = interval(end);
            tau = linspace(t0,tMax,Nt);

            % reachability cal.
            %[data, tau, extraOuts] = ...
            % HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
            [data, ~, ~] = HJIPDE_solve(-region.data, tau, schemeData,...
                'minVOverTime', HJIextraArgs); % careful with the reverse sign
            % now we flip the time and change the sign back
            dataCBF = -flip(data,4);
            data_new = squeeze(dataCBF(:,:,:,1));
            new_region.data = data_new;

        end

        if operatorSt(1) == 'G'
            new_region = region;
            interval = str2num(operatorSt(2:end));
            % reachability cal.
            tMax = interval(1);
            tau = linspace(t0,tMax,Nt);

            % reachability cal.
            %[data, tau, extraOuts] = ...
            % HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
            [data, ~, ~] = HJIPDE_solve(-region.data, tau, schemeData,...
                'minVOverTime', HJIextraArgs); % careful with the reverse sign
            % now we flip the time and change the sign back
            dataCBF = -flip(data,4);
            data_new = squeeze(dataCBF(:,:,:,1));
            new_region.data = data_new;
        end

        if operatorSt(1) == 'V'
            new_region = region;
            new_region.data = max(region.data,region2.data);
        end

        if operatorSt(1) == 'W'
            new_region = region;
            new_region.data = min(region.data,region2.data);
        end
    otherwise
        error('Type of region representation not supported.')
end

end

function new_region = setUnion(region,region2)
    distance = norm(region.c - region2.c);
    radius_err = abs(region.r - region2.r);
    radius_sum = region.r + region2.r;
    if distance<radius_err
        % almost concentric circles
        [~,ind] = max([region.r,region2.r]);
        if ind ==1 
            new_region = region;
        else
            new_region = region2;
        end  
    elseif distance>=radius_err && distance<=radius_sum
        % intersected, not fully contained
        warning('Overlapping. The new region is NOT a circular region. Treating as a circular one ANYWAY. ')
    else
        warning('Disjoint. The new region is NOT a circular region. Treating as a circular one ANYWAY. ')
    end

    [~,ind] = max([region.r,region2.r]);
    if ind ==1 
        new_region = region;
    else
        new_region = region2;
    end
end

function new_region = setIntersection(region,region2)
    distance = norm(region.c - region2.c);
    radius_err = abs(region.r - region2.r);
    radius_sum = region.r + region2.r;
    if distance<radius_err
        % almost concentric circles
        % choose smaller circle
        [~,ind] = max([region.r,region2.r]);
        if ind ==1 
            new_region = region2;
        else
            new_region = region;
        end  
    elseif distance>=radius_err && distance<=radius_sum
        % intersected, not fully contained
        warning('Overlapping. The new region is NOT a circular region. Treating as a circular one ANYWAY. ')
        largestIntersectingCircle = findLargestIntersectingCircle(region.c,region.r,region2.c,region2.r);
        new_region = region;
        new_region.r = largestIntersectingCircle.r;
        new_region.c = largestIntersectingCircle.c;
    else
        warning('Disjoint. The new region is NOT a circular region. Treating as a circular one ANYWAY. ')
        new_region = [];
    end

end

function largestIntersectingCircle = findLargestIntersectingCircle(center1, radius1, center2, radius2)
    % Calculate the distance between the two centers
    distance = norm(center1 - center2);
    
    % Check if the circles are disjoint
    if distance >= radius1 + radius2
        largestIntersectingCircle = [];
        return;
    end
    
    % Calculate the intersection points of the two circles
    d = distance;
    
    
    % Calculate the points of intersection
    p1 = center1 + radius1 * (center2 - center1) / d;
    p2 = center2 - radius2 * (center2 - center1) / d;
    
    
    % Calculate the radius of the largest intersecting circle
    largestRadius = norm(p1 - p2)/2;
    
    % Calculate the center of the largest intersecting circle
    largestCenter = (p1 + p2) / 2;
    
    % Store the results
    largestIntersectingCircle.c = largestCenter;
    largestIntersectingCircle.r = largestRadius;
end


