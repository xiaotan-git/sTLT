classdef operatorNodeObj < handle
    %OPERATORSETOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodeName % a string, like "Vee", "Wedge","F[a,b]","G[a,b]"
        nodeIndex
        id
        operatorNodeDuration
        startTimeInterval
    end
    
    methods
        function obj = operatorNodeObj(nodeName,nodeIndex)
            %OPERATORSETOBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.nodeName = nodeName;
            obj.nodeIndex = nodeIndex;
            
            % define operatorNodeDuration and startTimeInterval
            if strcmp(nodeName,'Vee') | strcmp(nodeName,'Wedge')
                obj.operatorNodeDuration = 0;
                obj.startTimeInterval = [0,0];
            end
            if strcmp(nodeName(1),'F')
                timeInterval = str2num(nodeName(2:end));
                obj.operatorNodeDuration = 0;
                obj.startTimeInterval = timeInterval;
            end
            if strcmp(nodeName(1),'G')
                timeInterval = str2num(nodeName(2:end));
                obj.operatorNodeDuration = timeInterval(2) - timeInterval(1);
                obj.startTimeInterval = [timeInterval(1) timeInterval(1)];
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
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

