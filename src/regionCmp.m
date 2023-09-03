function output = regionCmp(region1,region2,varargin)
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

% output: 0 - region1 = region2
%         1 - region1 contains region2
%         2 - region2 contains region1
%         3 - overlapping but no containment
%         4 - no overlapping


p = inputParser;
% % Define the expected parameters and their default values
% more details, see https://se.mathworks.com/help/matlab/ref/inputparser.html
% addParameter(p, 'paramName', defaultValue, validationFunction);
% optional propositional argument
% possible regionRep: circular, grid-data, levelset 
addOptional(p, 'regionRep', 'circular', @(x) ischar(x));
% optional name-value pair argument
% addParameter(p, 'umax', 1, @(x) isnumeric(x));
% addParameter(p, 'vRange', [-1 1], @(x) isnumeric(x));
% addParameter(p, 'wMax', 1, @(x) isnumeric(x));


parse(p, varargin{:});
% % parsedValue = p.Results.paramName;
% region = p.Results.region1;
% operatorSt = p.Results.operatorSt;
regionRep = p.Results.regionRep;
% umax = p.Results.umax;
% vRange = p.Results.vRange;
% wMax = p.Results.wMax;

switch regionRep
    case 'circular'
        distance = norm(region1.c - region2.c);
        radius_err = abs(region1.r - region2.r);
        % case 0
        if distance<10*eps && radius_err<10*eps
            output = 0;
            return
        end
        
        radius_sum = region1.r + region2.r;
        % case 1 or 2
        if distance<radius_err
            % almost concentric circles
            [~,ind] = max([region1.r,region2.r]);
            if ind ==1 
                output = 1;
            else
                output = 2;
            end
            return
        elseif distance>=radius_err && distance<=radius_sum
            % intersected, not fully contained
            output = 3;
            return
        else
            output = 4;
            return
        end

    case 'grid-data'
        % 
        if ~isequal(region1.grid,region2.grid)
            error('The grid is inconsistant for the two regions')
        end
        tf_1_contains_2 = all(region1.data(region2.data>=0)>=0);
        tf_2_contains_1 = all(region2.data(region1.data>=0)>=0);

        if tf_1_contains_2 && tf_2_contains_1
            output = 0;
            return
        else
            if tf_1_contains_2
                output = 1;
                return
            else
                if tf_2_contains_1
                    output = 2;
                    return
                else
                    % intersection
                    data = min(region1.data,region2.data);
                    tf_intersection = ~all(data<0);
                    if tf_intersection
                        output = 3;
                        return
                    else
                        output = 4;
                        return
                    end
                end
            end
        end
end




end
