function u_new = sat_func(u,lb,up)
%SAT_FUNC Summary of this function goes here
%   Detailed explanation goes here
if isempty(lb)
    lb = -Inf.*ones(size(u));
end

if isempty(up)
    up = Inf.*ones(size(u));
end

u_new = zeros(size(u));
for i=1:size(u,1)
    u_new(i) = min(max(lb(i),u(i)),up(i));
end
end
