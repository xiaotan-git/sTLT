function u = cbf_QP(A,b,lb,ub,unom)
%CBF_QP Summary of this function goes here
%   Detailed explanation goes here

%   solve    min\| u \|^2
%       s.t.  Au+b_vec >=0 %% following the convention in the paper

options = optimoptions('quadprog','Display','off');

switch nargin
    case 4
        m = size(A,2);
        H = diag(ones(m,1));
        f = zeros(m,1);
    case 5
        m = size(A,2);
        H = diag(ones(m,1));
        f = -unom;
end

[u,~,exitflag,~] = quadprog(H,f,-A,b,[],[],lb,ub,[],options);

if exitflag == -2
    error('The online QP is infeasible' )
end

end

