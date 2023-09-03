classdef setNodeObj < handle
    %STARTINGTIMEOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodeName % a string, like "X1"
        nodeIndex % index number
        id % a string, like 'X_0_1_2'
        setNodeDuration % a scalar
        startTimeInterval % a 2-D row vector, like [1,2]
        stlFormula % for example, root node -> full formula, predicate -> mu1 or an function
        region
    end
    
    methods
        function obj = setNodeObj(nodeName,nodeIndex)
            obj.nodeName = nodeName;
            obj.nodeIndex = nodeIndex;
        end
        
        function obj = setNodeObj_values(nodeName,startTimeInterval,setNodeDuration,region)
            %STARTINGTIMEOBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.nodeName = nodeName;
            obj.startTimeInterval = startTimeInterval;
            obj.setNodeDuration = setNodeDuration;
            obj.region = region;
        end
        
        function output = isInStartInterval(obj,t)
            if t >= obj.startTimeInterval(1) && t <= obj.startTimeInterval(2)
                output = true;
            else
                output = false;
            end
        end
        
%         function obj=updateStartInterval(obj,t,ts_bar_j)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             if inInterval(obj,t)
%                 obj.startTimeInterval = [t,t];
%             else
%                 time_shift = ts_bar_j - t;
%                 obj.startTimeInterval = obj.startTimeInterval- [time_shift,time_shift];
%             end
%         end
        
        function output = isInRegion(obj,x)
            switch obj.region.rep
                case 'circular'
                    if norm(x - obj.region.c) <=obj.region.r
                        output = true;
                    else
                        output = false;
                    end
                case 'grid-data'
                    val = eval_u(obj.region.grid,obj.region.data,x);
                    if val >= 0
                        output  = true;
                    else
                        output = false;
                    end
            end
        end
        
        function tf = eq(obj1,obj2)
            if obj1.nodeIndex == obj2.nodeIndex
                tf = true;
            else
                tf = false;
            end
        end
    end
end

