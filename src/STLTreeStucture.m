function [nodeList,parentList] = STLTreeStucture(phi)
%STLTREESTUCTURE Retrieve the sTLT tree structure from STL formula phi
% relevant properties: phi.type
%                         .id
%                         .phi
%                         .phi1
%                         .phi2
%                         .interval
%   Detailed explanation goes here

nodeList = {};
parentList = {};

phi_test_list = {phi};
setNodeCounter = 0;
while(~isempty(phi_test_list))
    % choose one formula to test and remove it from the test list 
    phi_test=  phi_test_list{1};
    phi_test_list(1)=[];
    
    % test 'or' 'and'
    if strcmp(phi_test.type,'or') || strcmp(phi_test.type,'and')
        nodeIndex = numel(nodeList)+1;
        nodeName = ['X_' num2str(setNodeCounter)];
        setNodeCounter = setNodeCounter+1;
        nodeId = ['X_' phi_test.id];
        setNodeTem = setNodeObj(nodeName,nodeIndex);
        setNodeTem.id = nodeId;
        setNodeTem.stlFormula = phi_test.st;
        
        nodeList{end+1} = setNodeTem;
        parentList{end+1} = -1;
        
        if strcmp(phi_test.type,'or')
            nodeType = 'Vee';
        else
            nodeType = 'Wedge';
        end
        
        nodeIndex = numel(nodeList)+1;
        nodeName = nodeType;
        nodeId = [nodeType(1) '_' phi_test.id];
        opNodeTem = operatorNodeObj(nodeName,nodeIndex);
        opNodeTem.id = nodeId;

        % set it now
        nodeList{end+1} = opNodeTem;
        parentList{end+1} = nodeIndex-1;
        
        % add two phis into the test list
        phi_test_list{end+1} = phi_test.phi1;
        phi_test_list{end+1} = phi_test.phi2;
    end
    
    if strcmp(phi_test.type,'eventually') || strcmp(phi_test.type,'always')
        nodeIndex = numel(nodeList)+1;
        nodeName = ['X_' num2str(setNodeCounter)];
        setNodeCounter = setNodeCounter+1;
        nodeId = ['X_' phi_test.id];
        setNodeTem = setNodeObj(nodeName,nodeIndex);
        setNodeTem.id = nodeId;
        setNodeTem.stlFormula = phi_test.st;
        
        nodeList{end+1} = setNodeTem;
        parentList{end+1} = -1;
        
        if strcmp(phi_test.type,'eventually')
            nodeType = 'F';
        else
            nodeType = 'G';
        end
        
        nodeName = [nodeType num2str(phi_test.interval)];
        nodeId = [nodeType '_' phi_test.id];
        nodeIndex = numel(nodeList)+1;
        opNodeTem = operatorNodeObj(nodeName,nodeIndex);
        opNodeTem.id = nodeId;
        
        nodeList{end+1} = opNodeTem;
        parentList{end+1} = nodeIndex-1;
        
        % add two phis into the test list
        phi_test_list{end+1} = phi_test.phi;
    end
    
    if strcmp(phi_test.type,'predicate')
        nodeIndex = numel(nodeList)+1;
        nodeName = ['X_' num2str(setNodeCounter)];
        setNodeCounter = setNodeCounter+1;
        nodeId = ['X_' phi_test.id];
        
        setNodeTem = setNodeObj(nodeName,nodeIndex);
        setNodeTem.id = nodeId;
        setNodeTem.stlFormula = phi_test.st;
        nodeList{end+1} = setNodeTem;
        % set later
        parentList{end+1} = -1;
    end
end

% now set the parent node for setNodes
for i = 1:numel(nodeList)
    if isa(nodeList{i},'setNodeObj')
        nodeIdi = nodeList{i}.id;

        % if it is the root node
        if strcmp(nodeIdi,['X_' phi.id])
            parentList{i} = 0;
            continue
        end

        % check for all other set nodes
        ind_nodeNamei = nodeIdi(2:end);
        
        for j = 1:numel(nodeList)
            nodeIdj = nodeList{j}.id;
            ind_nodeNamej = nodeIdj(2:end);
            if isa(nodeList{j},'operatorNodeObj') && ...
                    strcmp(ind_nodeNamei(1:end-2),ind_nodeNamej)
                parentList{i} = j;
            end
        end
    end
    

end
        


