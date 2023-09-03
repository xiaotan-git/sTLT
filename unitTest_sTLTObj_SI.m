% clc; clear all;
if case_numb ==1
    X0Node = setNodeObj('X',1);
    X1Node = setNodeObj('X1',2);
    X3Node = setNodeObj('X3',3);
    X5Node = setNodeObj('X5',4);
    X2Node = setNodeObj('X2',8);
    X4Node = setNodeObj('X4',9);
    X6Node = setNodeObj('X6',10);
    X7Node = setNodeObj('X7',11);
    X8Node = setNodeObj('X8',12);
    X9Node = setNodeObj('X9',13);

    Vee = operatorNodeObj('Vee',5);
    F_0_15 = operatorNodeObj('F[0,15]',6);
    G_2_10 = operatorNodeObj('G[2,10]',7);
    F1_0_15 = operatorNodeObj('F[0,15]',14);
    Wedge = operatorNodeObj('Wedge',15);
    G_0_10 = operatorNodeObj('G[0,10]',16);
    F_5_10 = operatorNodeObj('F[5,10]',17);

    nodeList = {X0Node,X1Node,X3Node,X5Node,Vee,F_0_15,G_2_10,X2Node,X4Node,X6Node,X7Node,X8Node,X9Node,F1_0_15,Wedge,G_0_10,F_5_10};
    parentNodeList = {[],Vee,F_0_15,G_2_10,X0Node,X1Node,X3Node,Vee,F1_0_15,Wedge,Wedge,G_0_10,F_5_10,X2Node,X4Node,X6Node,X7Node};
    % parentList1 = {0,5,6,7,1,2,3,5,14,15,15,16,17,8,9,10,11};
    parentList ={};
    for i = 1:length(parentNodeList)
        if ~isempty(parentNodeList{i})
            for j = 1:length(nodeList)
                if nodeList{j} == parentNodeList{i}
                    parentList{end+1} = j;
                end
            end
        else
            parentList{end+1} = 0;
        end
    end

    tree1 = sTLTObj(nodeList,parentList);
    tree1 = tree1.BUTran();
    tree1 = tree1.TDTran();
    % tree1.generate_completePath()
    tree1 = tree1.calStartTimeInterval();
    %tree1.updateStartTimeInterval(3,10);
    tree1.findAllTemporalFragments();

    draw_sTLT(tree1)

end

if case_numb == 2
    X0Node = setNodeObj('X',1);
    X1Node = setNodeObj('X1',2);
    X3Node = setNodeObj('X3',3);
    X5Node = setNodeObj('X5',4);
    X2Node = setNodeObj('X2',8);
    X4Node = setNodeObj('X4',9);
    X6Node = setNodeObj('X6',10);
    X7Node = setNodeObj('X7',11);
    X8Node = setNodeObj('X8',12);
    X9Node = setNodeObj('X9',13);

    Vee = operatorNodeObj('Vee',5);
    F_0_15 = operatorNodeObj('F[0,15]',6);
    G_2_10 = operatorNodeObj('G[2,10]',7);
    F1_0_15 = operatorNodeObj('F[0,15]',14);
    Wedge = operatorNodeObj('Wedge',15);
    G_0_10 = operatorNodeObj('G[0,10]',16);
    F_5_10 = operatorNodeObj('F[5,10]',17);

    nodeList = {X0Node,X1Node,X3Node,X5Node,Vee,F_0_15,G_2_10,X2Node,X4Node,X6Node,X7Node,X8Node,X9Node,F1_0_15,Wedge,G_0_10,F_5_10};
    parentNodeList = {[],Vee,F_0_15,G_2_10,X0Node,X1Node,X3Node,Vee,F1_0_15,Wedge,Wedge,G_0_10,F_5_10,X2Node,X4Node,X6Node,X7Node};
    % parentList1 = {0,5,6,7,1,2,3,5,14,15,15,16,17,8,9,10,11};
    parentList ={};
    for i = 1:length(parentNodeList)
        if ~isempty(parentNodeList{i})
            for j = 1:length(nodeList)
                if nodeList{j} == parentNodeList{i}
                    parentList{end+1} = j;
                end
            end
        else
            parentList{end+1} = 0;
        end
    end

    tree1 = sTLTObj(nodeList,parentList);
    tree1 = tree1.BUTran();
    tree1 = tree1.TDTran();
    % tree1.generate_completePath()
    tree1 = tree1.calStartTimeInterval();
    %tree1.updateStartTimeInterval(3,10);
    tree1.findAllTemporalFragments();

    draw_sTLT(tree1)
end
