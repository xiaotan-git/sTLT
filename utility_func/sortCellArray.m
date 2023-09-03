function sortedCellArray = sortCellArray(cellArray,sortDimen,sortOrder)
    % input cellArray = {[1,2],[3,5],[23,23],...}

    if nargin == 1
        sortDimen = 1;
        sortOrder = 'descend';
    else 
        if nargin == 2
            sortOrder = 'descend';
        end
    end

    arrayD = zeros(1,numel(cellArray));
    for i = 1:numel(cellArray)
        arrayD(i) =  cellArray{i}(sortDimen);
    end

    [~,id] = sort(arrayD,sortOrder);
    sortedCellArray = cellArray(:,id);
end