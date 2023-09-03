function [ind1,ind2] = isEleEqualCellArray(array1,array2)
%COMPARETWOCELLARRAY compare two cell arrays and return the index of elements with has a copy in another array
%   Ex. array1 = {1,2,3} = {[1]} {[2]} {[3]}; array2 = {2,3,4}; 
%  Function returns ind1= [0 1 1]; ind2 = [1 1 0];
if ~iscell(array1) || ~iscell(array2)
    % warning('Cell array expected for the function isEleEqualCellArray ')
    array1 = num2cell(array1);
    array2 = num2cell(array2);
end

ind1 = zeros(1,numel(array1)); ind2 = zeros(1,numel(array2)); 

for i = 1:numel(array1)
    for j = 1:numel(array2)
        if isequal(array1{i},array2{j})
            ind1(i) = 1;
            ind2(j) = 1;
        end
    end
end
end

