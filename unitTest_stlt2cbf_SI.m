switch case_numb
    case 1
        load('tree_singleIntegrator_example1.mat');
    case 2
        load('tree_singleIntegrator_example2.mat');
end


vMax = 1;
 
% cell array with double array as its entry
% {[ind_op, ind_set],[ind_op1, ind_set1],...}
branch_list = tree.findAllBranch();
if numel(branch_list)>1
    promt_str = ['Choose a branch out of ' num2str(numel(branch_list)) ': '];
    choose_branch = input(promt_str);
    if choose_branch>numel(tree.findAllBranch())
        error('Wrong branch number')
    end
    temporalFragments = tree.findAllTemporalFragments(choose_branch); 
else
    temporalFragments = tree.findAllTemporalFragments(); 
end


cbf_cell = {};timeInterval_cell = {};
% first we calculate timeIntervals [tb_low, ts_upp, tb_lower_child,tb_upp]
% this needs to be calculated from a button-up manner
temporalFragments = sortCellArray(temporalFragments,1,'descend');

for i = 1:numel(temporalFragments)
    operatorNode = tree.nodeList{temporalFragments{i}(1)};
    setNode = tree.nodeList{temporalFragments{i}(2)};
    % starting time interval
    ts_Interval = setNode.startTimeInterval;
    % Index of PA(PA(setNode))
    pa_setNodeIndex = tree.parentList{tree.parentList{temporalFragments{i}(2)}};
    pa_setNode = tree.nodeList{pa_setNodeIndex};
    ts_Interval_PA = pa_setNode.startTimeInterval;

    % time interval from Sec. III.D 2)
    % I prefer to use an earlier start time for the cbf
    tb_interval = [min(ts_Interval(1),ts_Interval_PA(2)+pa_setNode.setNodeDuration),...
        ts_Interval(2)+setNode.setNodeDuration ];
    % tb_interval = [min(ts_Interval(2),ts_Interval_PA(2)+pa_setNode.setNodeDuration),...
    %     ts_Interval(2)+setNode.setNodeDuration ];

    timeInterval = [tb_interval(1),ts_Interval(2), tb_interval(2)];
    c = setNode.region.c;
    r = setNode.region.r;
    obs = setNode.region.obstacle;

    % time interval for cbf builder [tb_low, ts_upp, tb_upp]
    % tb_low, tb_upp: lower and upper bound of the time domain of the barrier function
    % ts_upp: the upper bound time traj needs to be inside the set region
    % so for the time interval [tb_low, ts_upp], the system approaches the set region;
    % for the time interval [ts_upp, tb_upp], the system stays inside the set region

    % From condition 1 and 2 in Sec. III.D 2)
    cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
    cbf_i.temporalOperator = operatorNode.nodeName;
    cbf_i.setName = setNode.nodeName;
    cbf_cell{end+1} = cbf_i;

    use_childchildTemp = 0;
    if use_childchildTemp
        % find child temporal fragment
        childTempFrag = tree.findChildrenTemporalFragments(temporalFragments{i});
            
        if isempty(childTempFrag)
            % timeInterval = [ tb_lower(fi)  ts_upp(Xi)  tb_lower{f_children} tb_upp{fi}]
            % when at the leaf temporal fragments, set tb_lower{f_children}= tb_upp{fi}
            timeInterval = [tb_interval(1),ts_Interval(2), tb_interval(2)];
            c = setNode.region.c;
            r = setNode.region.r;
            obs = setNode.region.obstacle;

            % time interval for cbf builder [tb_low, ts_upp, tb_upp]
            % tb_low, tb_upp: lower and upper bound of the time domain of the barrier function
            % ts_upp: the upper bound time traj needs to be inside the set region
            % so for the time interval [tb_low, ts_upp], the system approaches the set region;
            % for the time interval [ts_upp, tb_upp], the system stays inside the set region

            % From condition 1 and 2 in Sec. III.D 2)
            cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
            cbf_i.temporalOperator = operatorNode.nodeName;
            cbf_i.setName = setNode.nodeName;
            cbf_cell{end+1} = cbf_i;
        else
            % get tb_lower{f_children} from childTempFrag
            tb_lower_child = inf; region_at_tb_lower_child = struct();
            for k = 1:numel(childTempFrag)
                nodeName_child = tree.nodeList{childTempFrag{k}{2}}.nodeName;
                for j = 1:numel(cbf_cell)
                    if strcmp(cbf_cell{j}.setName,nodeName_child) && ...
                        tb_lower_child>=cbf_cell{j}.timeInterval(1)
                        tb_lower_child = cbf_cell{j}.timeInterval(1);
                        region_at_tb_lower_child = cbf_cell{j}.regionAtTimeT(tb_lower_child);
                    end
                end
            end

            if isequal(tb_lower_child, tb_interval(2)) ||...
                (regionCmp(region_at_tb_lower_child,setNode.region,setNode.region.rep) ==0) ||...
                    (regionCmp(region_at_tb_lower_child,setNode.region,setNode.region.rep) ==1)
                % when the temporal fragment has no time overlap with its children temporal fragments
                % or when the region from the childen temporal fragments at tb_lower_child is equal to
                % or larger than the setNode region
                timeInterval = [tb_interval(1),ts_Interval(2),tb_interval(2)];
                timeInterval_cell{end+1} =  timeInterval;
                c = setNode.region.c;
                r = setNode.region.r;
                obs = setNode.region.obstacle;

                cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
                cbf_i.temporalOperator = operatorNode.nodeName;
                cbf_i.setName = setNode.nodeName;
                cbf_cell{end+1} = cbf_i;
            else
                % now tb_lower{f_children}< tb_upp{fi}
                % and the set node region is larger than its children's region_at_tb_lower_child
                timeInterval = [tb_interval(1),ts_Interval(2),tb_interval(2)];
                c = setNode.region.c;
                r = setNode.region.r;
                obs = setNode.region.obstacle;
                vMax = 1;
                cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
                cbf_i.temporalOperator = operatorNode.nodeName;
                cbf_i.setName = setNode.nodeName;
                cbf_cell{end+1} = cbf_i;

                % and we add one more cbf to enforce condition 3
                timeInterval = [tb_interval(1),tb_lower_child,tb_lower_child];
                c = region_at_tb_lower_child.c;
                r = region_at_tb_lower_child.r;
                obs = region_at_tb_lower_child.obstacle;
                cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
                cbf_i.temporalOperator = operatorNode.nodeName;
                cbf_i.setName = setNode.nodeName;
                cbf_cell{end+1} = cbf_i;

            end

        end
    end


end

% and then we check if some updates on the time interval is needed
% for i = 1:numel(temporalFragments)
%     % check if it has a child temporalFragment
%     % if yes, check if condition 3 holds; otherwise, no change
%     % if condition 3 does not hold, add one more barrier 
%     % with the time-interval [tb_low, tb_low_child, tb_low_child]
%     % and the final region is the starting region of the child temporal
%     % fragment
%     
%     timeInterval = timeInterval_cell{i};
%     c = setNode.region.c;
%     r = setNode.region.r;
%     obs = setNode.region.obstacle;
%     vMax = 1;
%     cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
%     cbf_i.temporalOperator = operatorNode.nodeName;
%     cbf_i.setName = setNode.nodeName;
%     
%     cbf_cell{end+1} = cbf_i;
% end
switch case_numb
    case 1
        save('cbf_cell_singleIntegrator_example1.mat','cbf_cell');
    case 2
        save('cbf_cell_singleIntegrator_example2.mat','cbf_cell');
end

