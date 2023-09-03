clc; clear all;close all;

% st5 = 'G[0,1] F[2,3] mu1 Vee alw_[0,10] (x>=3)';
% phi1 = STLformula('0',st5);
% [nodeList,parentList] = STLTreeStucture(phi1);
% tree1 = sTLTObj(nodeList,parentList);
% tree1.draw_sTLT()

case_numb = 0;
sys_dyn = 'unicycle';
vRange = [-1 1];
wMax = 1;
use_valueFunction = true;

% inputs
switch case_numb
    case 0
        % trivial case STL
        st = 'F[0,15] mu1';
        predicateList = {'mu1','mu2', 'mu3'}';
        R1.c = [-4;-4]; R1.r = 1; R1.obstacle = false;
        
        % Create a containers.Map
        predMap = containers.Map;
        predMap('mu1') = R1;

        if use_valueFunction
            grid_min = [-8;-8;-pi];
            grid_max = [8;5;pi];
        end
        
    case 1
        % STL
        st = 'F[0,15] G[2,10] mu1 Vee F[0,15](G[0,10] mu2 Wedge F[5,10] mu3 )';
        predicateList = {'mu1','mu2', 'mu3'}';
        R1.c = [-4;-4]; R1.r = 1; R1.obstacle = false;
        R2.c = [4;0];R2.r = 4;R2.obstacle = false;
        R3.c = [1;-4];R3.r = 2;R3.obstacle = false;

        % Create a containers.Map
        predMap = containers.Map;
        predMap('mu1') = R1;
        predMap('mu2') = R2;
        predMap('mu3') = R3;

        if use_valueFunction
            grid_min = [-8;-8;-pi];
            grid_max = [8;5;pi];
        end
       
% case 2
    case 2
        st = ['G[0,1] F[2,3] mu1 and F[6,7] G[1,2] mu2' ...
        ' and F[15,16](G[0,4] mu3 and F[1,4] mu1) and G[0,20] (not mu4)'];
        predicateList = {'mu1', 'mu2', 'mu3', 'mu4'}';
        R1.c = [-2;2]; R1.r = 1; R1.obstacle = false;
        R2.c = [1.5;0.5];R2.r = 1;R2.obstacle = false;
        R3.c = [-1;0];R3.r = 1.5;R3.obstacle = false;
        R4.c = [0;1.5]; R4.r = 1;R4.obstacle = false;
        % Create a containers.Map
        predMap = containers.Map;
        predMap('mu1') = R1;
        predMap('mu2') = R2;
        predMap('mu3') = R3;
        predMap('mu4') = R4;
        
        if use_valueFunction
            grid_min = [-5; -5; -pi]; % Lower corner of computation domain
            grid_max = [5; 5; pi];    % Upper corner of computation domain
        end
end

%% offline computation
% parse the string formula to an STL formula obj
phi = STLformula('0',st);
% construct the tree structure
[nodeList,parentList] = STLTreeStucture(phi);
% construct the sTLT tree
tree = sTLTObj(nodeList,parentList);
% illustrate the tree
tree.draw_sTLT()

tree.BUTran(); % botton up transverse
tree.TDTran(); % top-down transverse

% assigning relevant setNodeObj its region/value function
if use_valueFunction
    predValueMap = containers.Map;
    keyArray = keys(predMap);
    for i = 1:numel(keyArray)
        % [grid,data] = build_unicycleCBF_helper_cylinder(grid_min,grid_max,center,radius)
        [grid,data] = build_unicycleCBF_helper_cylinder(grid_min,grid_max,predMap(keyArray{i}).c,predMap(keyArray{i}).r);
        region.grid = grid; region.data = data;region.rep = 'grid-data';
        predValueMap(keyArray{i}) = region;
    end
end

if use_valueFunction
    tree.assignPredicateRegion(predValueMap);
end

tree.updateALLRegions();
% 
% % Algorithm 3
tree.calStartTimeInterval();
% % Algorithm 4 
% tree1.updateStartTimeInterval(nodeIndex,updateTime)
% tree.updateStartTimeInterval(9,10); 
% % for cbf construction
tree.findAllTemporalFragments();
tree.isAllVeeTopLayer();
tree.findAllBranch();

switch case_numb
    case 1
        save('tree_unicycle_example1.mat','tree')
    case 2
        save('tree_unicylce_example2.mat','tree')
end


