function [grid,data] = build_unicycleCBF_helper_rectangle(grid_min,grid_max,center,radius)
%BUILD_UNICYCLECBF_HELPER Summary of this function goes here
%   Detailed explanation goes here


% Grid
N = [41; 41; 41];         % Number of grid points per dimension
pdDims = 3;               % 3rd dimension is periodic
grid = createGrid(grid_min, grid_max, N, pdDims);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions


% target set
% data = shapeRectangleByCenter(grid, center, widths)
% CBF construction for the set X = {x: data0(x)>=0} 
center = [center;0]; % 3d
widths = [radius;radius;inf]; % 3d
data0 = -shapeRectangleByCenter(grid, center, widths);

data = data0;
end

