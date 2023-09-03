function plot_traj_color(x_vec,t_vec)
%PLOT_TRAJ_COLOR Summary of this function goes here
%   Detailed explanation goes here
hold on;
loops = 100;
sample = ceil(size(x_vec,2)/loops);
x = x_vec(1,1:sample:end); y = x_vec(2,1:sample:end); z = zeros(size(x));
t = t_vec(1,1:sample:end);

% plot(x_vec(1,1),x_vec(2,1),'ko','MarkerSize',12);

col = t;  % This is the color, vary with t in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    
end

