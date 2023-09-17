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


cbf_cell = {};
for i = 1:numel(temporalFragments)
    operatorNode = tree.nodeList{temporalFragments{i}(1)};
    setNode = tree.nodeList{temporalFragments{i}(2)};
    ts_Interval = setNode.startTimeInterval;
    % run algorithm 3 if not done before
    if isempty(ts_Interval)
        tree.calStartTimeInterval();
    end
    pa_setNodeIndex = tree.parentList{tree.parentList{temporalFragments{i}(2)}};
    pa_setNode = tree.nodeList{pa_setNodeIndex};
    ts_Interval_PA = pa_setNode.startTimeInterval;

    % time interval from Sec. III.D 2)
    tb_interval = [min(ts_Interval(1),ts_Interval_PA(2)+pa_setNode.setNodeDuration),...
        ts_Interval(2)+setNode.setNodeDuration ];

    % time interval for cbf builder [tb_low, tset_upp, tb_upp]
    % Based on CBF design principle 2, the state needs to stay in the set
    % node region during [tset_upp, tb_upp]    
    timeInterval = [tb_interval(1),ts_Interval(2),tb_interval(2)];
    c = setNode.region.c;
    r = setNode.region.r;
    obs = setNode.region.obstacle;
    cbf_i = singleIntegratorCBF(timeInterval,c,r,vMax,obs);
    cbf_i.temporalOperator = operatorNode.nodeName;
    cbf_i.setName = setNode.nodeName;
    cbf_cell{end+1} = cbf_i;
end
switch case_numb
    case 1
        save('cbf_cell_singleIntegrator_example1.mat','cbf_cell');
    case 2
        save('cbf_cell_singleIntegrator_example2.mat','cbf_cell');
end

