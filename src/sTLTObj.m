classdef sTLTObj < handle
    %STLTOBJ Summary of this class goes here
    %   Detailed explanation goes here
    % by list, I mean cell array
    
    properties
        nodeList %{X0 X1 X2 X3 X4 .. Vee .. Wedge F[a,b] ... G[a,b] }
        nodeIndexList % [1,2,3,...,M]
        setNodeList %{X0 X1 X3 ... XN}%}
        setNodeIndex
        operatorNodeList
        operatorNodeIndex
        parentList %{0, 3, 5, ...} 0 -- No parent set
        childrenList %{[1,2], [2],...}
        topDownTransverseIndex % index list for node i to do a top-down transverse
        bottomUpTransverseIndex % index list for node i to do a bottom-up transverse
        temporalFragments % {[5,2],[10,13],...} where [5,2] is one temporal fragment
    end
    
    methods
        function obj = sTLTObj(nodeList,parentList)
            %STLTOBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.nodeList = nodeList;
            obj.nodeIndexList = [1:length(nodeList)];
            setNodeList = {};
            setNodeIndex = [];
            operatorNodeList = {};
            operatorNodeIndex = [];
            for i = 1:length(nodeList)
                if isa(nodeList{i},'setNodeObj')
                    setNodeList{end+1}= nodeList{i};
                    setNodeIndex(end+1) = i;
                end
                if isa(nodeList{i},'operatorNodeObj')
                    operatorNodeList{end+1}= nodeList{i};
                    operatorNodeIndex(end+1) = i;
                end
            end
            
            obj.setNodeList=setNodeList;
            obj.setNodeIndex = setNodeIndex;
            obj.operatorNodeList = operatorNodeList;
            obj.operatorNodeIndex = operatorNodeIndex;
            
            obj.parentList = parentList;
            obj.childrenList = {};
            for i = 1:length(nodeList)
                i_childrenList = [];
                for j = 1:length(parentList)
                   if i == parentList{j}
                       i_childrenList(end+1)=j;
                   end
                end
                obj.childrenList{end+1} = i_childrenList;
            end            
            
        end
        
        function output = findAllParent(obj,i)
            % gives all parent nodes of node i in the sTLT
            if obj.parentList{i}~=0
                parent_node = obj.parentList{i};
                output = [findAllParent(obj,parent_node), parent_node];
            else
                output = [];
            end
        end
        
        function output = findAllChildren(obj,nodes2Check)
            % gives all child nodes of node i in the sTLT
            flag = 0; childrenNodes=[];
            for i = 1:length(nodes2Check)
                nodeIndex = nodes2Check(i);
                if ~isempty(obj.childrenList{nodeIndex})
                    children_node_i = obj.childrenList{nodeIndex};
                    childrenNodes = [childrenNodes children_node_i];
                    flag = 1;
                end
            end
            if flag == 1
                output = [childrenNodes findAllChildren(obj,childrenNodes)];
            else
                output = [];
            end
        end
        
        
        function obj = BUTran(obj)
            bottomUpTransverseIndex = {};
            for i = 1: length(obj.nodeList)
                i_BUTranIndex = findAllParent(obj,i);
                bottomUpTransverseIndex{end+1} = i_BUTranIndex;
            end
            obj.bottomUpTransverseIndex = bottomUpTransverseIndex;
        end
        
        function obj = TDTran(obj)
            topDownTransverseIndex={};
            for i = 1:length(obj.nodeList)
                i_TDTranIndex = findAllChildren(obj,i);
                topDownTransverseIndex{end+1} = i_TDTranIndex;
            end
            obj.topDownTransverseIndex = topDownTransverseIndex;
        end
        
        function temporalFragments = findAllTemporalFragments(obj,varargin)
            % Not finished
            temporalFragments = {};
            switch numel(varargin)
                case 0 % default case
                    node2check = obj.nodeIndexList;
                case 1 % findAllTemporalFragments(obj,complete_path)
                    if  isnumeric(varargin{1})
                        choose_branch  = varargin{1};
                        branch_list = findAllBranch(obj);
                        complete_paths = branch_list{choose_branch}; % possible several paths
                        node2check = [];
                        for i = 1:numel(complete_paths)
                            [ind1,~] = isEleEqualCellArray(complete_paths{i},node2check);
                            node2add = complete_paths{i}(ind1 ==0);
                            node2check = [node2check,node2add];
                            % disp(node2check)
                        end
                    else
                        node2check = varargin{1};
                    end                
            end
            % disp(node2check)
            for i = 1:length(node2check)
                node_i_ind = node2check(i);
                node_i = obj.nodeList{node_i_ind};
                if isa(node_i,'setNodeObj')
                    parentNodeIndex = obj.parentList{node_i_ind};
                    if parentNodeIndex ~=0
                        parentNode = obj.nodeList{parentNodeIndex};
                        parentNodeName = parentNode.nodeName;
                        if strcmp(parentNodeName(1),'F') || strcmp(parentNodeName(1),'G')
                            temporalFragments{end+1} = [parentNodeIndex,node_i_ind];
                        end
                    end
                end
            end

            obj.temporalFragments = temporalFragments;    
        end

        function subtree = subtree(obj,nodeIndex)
            if isa(obj.nodeList{nodeIndex},'operatorNodeObj')
                nodeIndex = obj.parentList{nodeIndex}; % nodeIndex is always a set node
            end

            subtree_nodeList = obj.findAllChildren(nodeIndex);
            subtree_nodeList = [nodeIndex,subtree_nodeList]; % add itself
            subtree_node = {};
            for i = 1:numel(subtree_nodeList)
                subtree_node{end+1} = obj.nodeList{subtree_nodeList(i)};
            end
            subtree_parentList = {0};
            for i = 2:numel(subtree_nodeList)
                for j = 1:numel(subtree_nodeList)
                    if obj.parentList{subtree_nodeList(i)} == subtree_nodeList(j)
                        subtree_parentList{end+1} = j;
                    end
                end
            end
            subtree = sTLTObj(subtree_node,subtree_parentList);

        end

        function topTemFrags = findAllTopTemporalFragments(obj)
            % return all top-level temporal fragments, i.e., without predecessors
            allTempFrags = obj.findAllTemporalFragments();
            topTemFrags = {};
            for i = 1:numel(allTempFrags)
                % if there is no predecessor
                temporalNodeIndex = allTempFrags{i}(1);
                parentNodes = obj.findAllParent(temporalNodeIndex);
                toplevel = 1;
                for j = 1:numel(parentNodes)
                    if isa(obj.nodeList{parentNodes(j)},'operatorNodeObj')
                        if regexp(obj.nodeList{parentNodes(j)}.nodeName,'[GF]')
                            toplevel = 0;
                        end
                    end
                end
                if toplevel
                    temFrag_i = allTempFrags{i};
                    topTemFrags{end+1}= temFrag_i;
                end
            end

        end
        
        function children_temFrag = findChildrenTemporalFragments(obj,temFrag)
            % return the immediate child temporal fragments of a given one
            % input: [8,9] as the index in nodeList
            % output: {[10,11], [12,13],...} or {} if no child temporal fragment exists
            subSTLformula = obj.nodeList{temFrag(2)}.stlFormula;
            % need to fix STL2sTLT
            if iscell(subSTLformula)
                subSTLformula = subSTLformula{:};
            end
            if all(regexp(subSTLformula, '(alw|ev)'))
                childTempFrag_Flag = 1;
            else
                childTempFrag_Flag = 0;
            end

            if childTempFrag_Flag
                children_temFrag = {};
                subtree = obj.subtree(temFrag(2));
                topLayerTemFrags_subtree = subtree.findAllTopTemporalFragments();
                % now retreive it back, matching through set node name
                for i = 1:numel(topLayerTemFrags_subtree)
                    setNodeName = subtree.nodeList{topLayerTemFrags_subtree{i}(2)}.nodeName;
                    for j = 1:numel(obj.setNodeIndex)
                        if strcmp(setNodeName,obj.nodeList{obj.setNodeIndex(j)}.nodeName)
                            setNodeIndex_j = obj.setNodeIndex(j);
                            temFrag_j = [obj.parentList(setNodeIndex_j), setNodeIndex_j];
                            children_temFrag{end+1}= temFrag_j;
                        end
                    end
                end

            else
                children_temFrag = {};
            end

        end


        function draw_sTLT(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            matlabTree = convert2MatlabDigragh(obj);
            figure
            p = plot(matlabTree,'Marker',matlabTree.Nodes.nodeMarker,'MarkerSize',30,'NodeColor',[0.90,0.90,0.90],'Layout','layered');
                % 'NodeFontSize',14,'linewidth', 1,...
                % 'NodeColor',[0.90,0.90,0.90],'NodeFontWeight','bold',...
            p.NodeLabel = {};   
            text(p.XData,p.YData,matlabTree.Nodes.nodeLabel, ...  %the first 2 arguments provide the position and we can give individual X and Ypositions of each node
            'VerticalAlignment','middle',...                
            'HorizontalAlignment', 'center',... % alignment with the node
            'FontSize', 12)
            box off;
            axis off;
        end
        
        function matlabTree = convert2MatlabDigragh(obj)
            adjMatrix = zeros(length(obj.nodeList),length(obj.nodeList));
            for i = 1:length(obj.nodeList)
                for j = 1:length(obj.nodeList)
                    if j == obj.parentList{i}
                        adjMatrix(j,i) = 1;
                    end
                end
            end
            matlabTree = digraph(adjMatrix);
            
            nodeLabel = {}; nodeMarker = {};
            for i = 1:length(obj.nodeList)
                if regexp(obj.nodeList{i}.nodeName,'\<X_')
                    label = ['X_{' obj.nodeList{i}.nodeName(3:end) '}' ];
                    nodeLabel{end+1} = label;
                else 
                    if regexp(obj.nodeList{i}.nodeName,'Vee|Wedge')
                        label = ['\' lower(obj.nodeList{i}.nodeName)];
                        nodeLabel{end+1} = label;
                    else
                        nodeLabel{end+1} = obj.nodeList{i}.nodeName;
                    end
                end
                if isa(obj.nodeList{i},'setNodeObj')
                    nodeMarker{end+1} = 'square';
                else
                    nodeMarker{end+1} = 'o';
                end
            end 
            nodeLabel = nodeLabel'; nodeMarker = nodeMarker';
            matlabTree.Nodes.nodeLabel = nodeLabel;
            matlabTree.Nodes.nodeMarker = nodeMarker;
        end
        
        function completePathList = generate_completePath(obj)
            leafNodeList = {};
            for i=1:length(obj.nodeList)
                if isempty(obj.childrenList{i})
                    leafNodeList{end+1} = i;
                end
            end
            completePathList = {};
            for i = 1:length(leafNodeList)
                completePathList{end+1} = [findAllParent(obj,leafNodeList{i}) leafNodeList{i}];
            end
        end
        
%         function 
        
        function obj = calStartTimeInterval(obj) % Alg 3
            % cal obj.setNodeList{i}.startTimeInterval;
            for i=1:length(obj.nodeList)
                if obj.parentList{i}==0
                    RootIndex = i;
                    obj.nodeList{i}.startTimeInterval = [0, 0];
                    obj.nodeList{i}.setNodeDuration = 0;
                end
            end
            
            for j = 1:length(obj.topDownTransverseIndex{RootIndex})
                nodeIndex = obj.topDownTransverseIndex{RootIndex}(j);
                node_j = obj.nodeList{nodeIndex};
%                 disp(class(node_j))
                
                if isa(node_j,'setNodeObj') %setNodeIndex
                    parentNodeIndex = obj.parentList{nodeIndex};
                    parentParentNodeIndex = obj.parentList{parentNodeIndex};
                    obj.nodeList{nodeIndex}.startTimeInterval = ...
                        obj.nodeList{parentParentNodeIndex}.startTimeInterval ...
                        + obj.nodeList{parentNodeIndex}.startTimeInterval;
                    obj.nodeList{nodeIndex}.setNodeDuration = ...
                        obj.nodeList{parentParentNodeIndex}.setNodeDuration + ...
                        + obj.nodeList{parentNodeIndex}.operatorNodeDuration;
                end
            end
            
        end
        
        function set2UpdateIndexList = updateStartTimeInterval(obj,relativeRootNodeIndex,t) % Alg. 4
            if relativeRootNodeIndex>numel(obj.nodeList) || ~isa(obj.nodeList{relativeRootNodeIndex},'setNodeObj')
                error('Check the node index before updating start time interval')
            end
            obj.nodeList{relativeRootNodeIndex}.startTimeInterval = [t, t];
            
            set2UpdateIndexList= {relativeRootNodeIndex};
            for j = 1:length(obj.topDownTransverseIndex{relativeRootNodeIndex})
                nodeIndex = obj.topDownTransverseIndex{relativeRootNodeIndex}(j);
                node_j = obj.nodeList{nodeIndex};
%                 disp(class(node_j))
                
                if isa(node_j,'setNodeObj') %setNodeIndex
                    parentNodeIndex = obj.parentList{nodeIndex};
                    parentParentNodeIndex = obj.parentList{parentNodeIndex};
                    obj.nodeList{nodeIndex}.startTimeInterval = ...
                        obj.nodeList{parentParentNodeIndex}.startTimeInterval ...
                        + obj.nodeList{parentNodeIndex}.startTimeInterval;
                    
                    set2UpdateIndexList{end+1} = nodeIndex;
                end
            end
        end
        
        function setIndex = findSetNode(obj,specificSetName)
            setIndex = 0;
            for i = 1:length(obj.nodeList)
                if strcmp(obj.nodeList{i}.nodeName,specificSetName)
                    setIndex = i;
                end
            end
        end
        
        function leafNodes = findLeafNode(obj)
            % giving out an array [1,2,3,...]
            leafNodes = [];
            for i = 1:length(obj.nodeList)
                if isempty(findAllChildren(obj,i))
                    leafNodes(end+1) = i;
                end
            end
        end
        
        function obj = assignPredicateRegion(obj,predMap)
            % this function is to map the predicate-region pair defined in a
            % containers.Map to the corresponding setNode
            leafNodes = findLeafNode(obj);
            for i = 1:numel(leafNodes)
                st = obj.nodeList{leafNodes(i)}.stlFormula;
                pred_st = regexprep(st, '\s+|not|\(|\)', '');
                region = predMap(pred_st);
                % check for obstacle region
                % here we treat the cases 'mu1','not mu1','not not mu1'
                % since we assumed positive normal form (PNF)
                if mod(numel(regexp(st,'\<not\>','match')),2)
                    % when it is an obstacle region like 'not mu1', 'not not not mu1'
                    region.obstacle = ~region.obstacle;
                    if strcmp(region.rep,'grid-data')
                        region.data = - region.data;
                    end
                end
                obj.nodeList{leafNodes(i)}.region = region;
            end
            
        end
        
        function obj = updateALLRegions(obj,varargin)
            % this function is to calculate the regions in all temporal
            % fragments in the sTLT tree. It is assumed that the regions of
            % all the leaf nodes are well-defined already
            switch numel(varargin)
                case 0
                    leafNodes = findLeafNode(obj); % cell array
                    setNodeIndex = obj.setNodeIndex; % 
                    setNode2Update = setdiff(setNodeIndex,leafNodes);
                case 1
                    leafNodes = varargin{1};
                    setNodeIndex = obj.setNodeIndex;
                    setNode2Update = setdiff(setNodeIndex,leafNodes);
                case 3
                    leafNodes = varargin{1};
                    setNodeIndex = varargin{2};
                    setNode2Update = setdiff(setNodeIndex,leafNodes);
            end
            
            for i = 1:numel(leafNodes)
                buTran = obj.bottomUpTransverseIndex{leafNodes(i)};  %[1,2,3,...]
                buTran = flip(buTran); % So it starts with the nearest parent
                
                for j = 1:numel(buTran)
                    node_j = obj.nodeList{buTran(j)};
                    % find an operator node and its child setnode
                    if isa(node_j,'operatorNodeObj') %operatorNode Index
                        pa_j = buTran(j+1);
                        if j >1
                            ch_j = buTran(j-1);
                            ch_node = obj.nodeList{ch_j};
                        else
                            ch_node = obj.nodeList{leafNodes(i)};
                        end
                        
                        operator = node_j.nodeName;
                        region = ch_node.region;
                        % Now calculate the backward reachable set
                        if operator(1)=='F' || operator(1)=='G'
                            pa_node = obj.nodeList{pa_j};
                            % disp('the index of the set node to be updated F or G: ')
                            % disp(pa_j)
                            if ~isempty(region)
                                pa_node.region = regionCal(region,operator,...
                                'regionRep',region.rep);
                                % remove it from setNode2Update
                                setNode2Update(setNode2Update == pa_j) = [];
                            end
                        else
                            ch_list = obj.childrenList{buTran(j)};
                            ch1_node = obj.nodeList{ch_list(1)};
                            ch2_node = obj.nodeList{ch_list(2)};
                            if ~isempty(ch1_node.region) && ~isempty(ch2_node.region)
                                region1 = ch1_node.region;
                                region2 = ch2_node.region;
                                pa_node = obj.nodeList{pa_j};
                                % disp('the index of the set node to be updated V or W: ')
                                % disp(pa_j)
                                pa_node.region = regionCal(region1,operator,region2,...
                                'regionRep',region.rep);
                                % remove it from setNode2Update
                                setNode2Update(setNode2Update == pa_j) = [];
                            end
                        end
                    end
                end
                              
                                    
            end

            if ~isempty(setNode2Update)
                new_leafnode = {};
                setNodeUpdated = setdiff(obj.setNodeIndex,setNode2Update);
                for i = 1:numel(setNodeUpdated)
                    setNodei_Index = setNodeUpdated(i);
                    PAPA_setNodei_Index = node.parentList{node.parentList{setNodei_Index}};
                    if ~isempty(obj.nodeList{setNodei_Index}.region) && ...
                        isempty(obj.nodeList{PAPA_setNodei_Index}.region)
                        new_leafnode{end+1} = setNodei_Index;
                    end
                end

                obj = updateALLRegions(new_leafnode,setNode2Update);
            end

        end

        function tf = isAllVeeTopLayer(obj)
            % return true or false on whether vees only in top layers
            % When no vee operator exists, return true
            tf = true;
            for i = 1:numel(obj.operatorNodeIndex)
                nodeInd = obj.operatorNodeIndex(i);
                if strcmp(obj.nodeList{nodeInd}.nodeName,'Vee')
                    parent_i_List = obj.findAllParent(nodeInd); % an array
                    for j = 1:numel(parent_i_List)
                        node_j_Ind = obj.operatorNodeIndex(parent_i_List(j));
                        node_j = obj.nodeList{node_j_Ind};
                        if isa(node_j,'operatorNodeObj') && (~strcmp(node_j.nodeName,'Vee'))
                            tf = false;
                            return;
                        end

                    end
                end
            end
        end

        function branchList = findAllBranch(obj)
            % return a cell array with each element consisting of one or several complete pathes
            % e.g. branchList{1} = [1,2,3...]
            branchList = {};
            if ~isAllVeeTopLayer(obj)
                warning('The sTLT does NOT contain all Vee operators in its top layers. Check again.');
                draw_sTLT(obj);
                return;
            end
            
            % check if Vee exists
            % if no Vee exists, then all the complete paths are the same group
            num_vee = 0;
            for i = 1:numel(obj.operatorNodeList)
                node_i = obj.operatorNodeList{i};
                if strcmp(node_i.nodeName,'Vee')
                    num_vee = num_vee+1;
                end
            end
            if num_vee==0
                completePathList = generate_completePath(obj);
                branchList = {completePathList};
                return;
            end

            % now Vee operator exists and lies in top layer
            completePathList = generate_completePath(obj);
            % construct first branch that consists of the first path
            path_l = completePathList{1}; % letter l not number 1
            branch_l = {path_l}; vee_ind = []; k_l = 0;
            for i = 1:numel(path_l)
                node_i = obj.nodeList{path_l(i)};
                if isa(node_i,'operatorNodeObj') && strcmp(node_i.nodeName,'Vee') && k_l <= i
                    k_l = i;
                end
            end
            branchList{end+1} = branch_l;
            vee_ind(end+1) = k_l;

            % for other paths, check if it belongs to any existing branches. Otherwise, creat a new branch
            for f = 2:numel(completePathList)
                path_f = completePathList{f};
                k_f = 0;
                for i = 1:numel(path_f)
                    node_i = obj.nodeList{path_f(i)};
                    if isa(node_i,'operatorNodeObj') && strcmp(node_i.nodeName,'Vee') && k_f <= i
                        k_f = i;
                    end
                end

                % check if it matches existing branches
                % if yes, add it to the existing one
                % otherwise, create a new branch
                match_flag = 0;
                for i = 1:numel(branchList)
                    branch_i = branchList{i};
                    cpl_path_i = branch_i{1}; % take the first complete path in branch i
                    % disp(cpl_path_i)
                    % disp(vee_ind)
                    if (vee_ind(i) == k_f) && isequal(path_f(1,1:k_f+1),cpl_path_i(1,1:k_f+1))
                        % noting here we also check for the set node X_k_f
                        branchList{i}{end+1} = path_f;
                        match_flag = 1;
                    end
                end
                
                if match_flag == 0
                    branchList{end+1} =  {path_f};
                    vee_ind(end+1) = k_f;
                end

            end
            

        end
        
    end
end

