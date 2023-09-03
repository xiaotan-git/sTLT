function plot_unicycle_oriented(x,gain)
%PLOT_UNICYCLE_ORIENTED Summary of this function goes here
%   Detailed explanation goes here
position = x(1:2,1);
orientation = x(3,1);

% gain = 1.0;

x1_0 = gain*[0.2;0]; x2_0 = gain*[0;0.1]; x3_0 = gain*[0;-0.1];
x4_0 = gain*[-0.2;0.1];x5_0 = gain*[-0.2;-0.1];
R = [cos(orientation) -sin(orientation); sin(orientation) cos(orientation)];

x1 = position + R*x1_0;
x2 = position + R*x2_0;
x3 = position + R*x3_0;
x4 = position + R*x4_0;
x5 = position + R*x5_0;

hold on;
plot([x1(1) x2(1)],[x1(2) x2(2)],'k');
plot([x3(1) x2(1)],[x3(2) x2(2)],'k');
plot([x1(1) x3(1)],[x1(2) x3(2)],'k');
plot([x2(1) x4(1)],[x2(2) x4(2)],'k');
plot([x5(1) x3(1)],[x5(2) x3(2)],'k');
plot([x4(1) x5(1)],[x4(2) x5(2)],'k');

% plot(position(1),position(2),'o');
axis equal
end

