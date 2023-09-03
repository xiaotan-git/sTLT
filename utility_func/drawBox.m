function drawBox(centers,radii,fillFlag)
%DRAWCIRCLE Summary of this function goes here
%   Detailed explanation goes here
if nargin ==2
    fillFlag = 0;
end

% Q = radius^2.*diag(ones(size(center)));
% q = center;
% E = ellipsoid(Q,q);
% linObj = plot(E,[1,2]);
% if fillFlag == 1
%     linObj.Color = 'none';
%     p = patch(linObj.XData,linObj.YData,'red');
%     p.FaceColor =  [0.6,0.2,0];
%     p.EdgeColor = 'none';
%     p.FaceAlpha = 0.5;
% end
if size(centers,2)~=2
    centers = centers';
end
if size(radii,2)~=1
    radii = radii';
end

for i = 1:size(centers,1)
    rec = rectangle('Position', [centers(1, 1) - radii, centers(1, 2) - radii, 2*radii, 2*radii], ...
            'EdgeColor', 'b', 'LineWidth', 1);
        
    if fillFlag ==1
        rec.EdgeColor = 'None';
        rec.FaceColor = [0.8,0.8,0.8];
    end
end

end

