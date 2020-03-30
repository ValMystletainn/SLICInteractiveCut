function newL = postProcess(L)
% post process of a label image
% it will enhanced the connectivity of the input L. Results in the single
% connected area of each label
% Input:
% L -- the raw label image
% Output:
% newL -- enhanced the connectivity patch of input - L
    num = max(L(:));
    newL = -1 * ones(size(L));
    [M, N] = size(L);
    % cut out the area that not the biggest connected area in each labels
    for i = 1 : num
        mask = L == i;
        connectedGroup = bwlabel(mask);
        stats = regionprops(connectedGroup);
        areas = cat(1, stats.Area);
        ind = find(areas ==max(areas));
        ind = ind(1);
        mask(find(connectedGroup~=ind))=0;
        newL(mask) = i;
    end
    % fill the hole of each label to be a Compact set
    for i = 1 : num
            mask = newL == i;   
            mask = imfill(mask,'holes'); 
            newL(mask) = i;
    end
    
    % assign those -1 point by the nearest, not -1, left/top pixel
    indNoconnect = find(newL == -1);
    for i = 1 : length(indNoconnect)
        noConnectedInd = indNoconnect(i);
        if (noConnectedInd == 1)
            newL(noConnectedInd) = 1;
            continue
        end
        [m, n] = ind2sub(size(L), noConnectedInd);
        searchsize = 2;
        breakpoint = max(m, n);
        while(searchsize < breakpoint)
            top = max(1, m - searchsize + 1);
            left = max(1, n - searchsize + 1);
            % create the search index
            searchArea = [fliplr(newL(top, left:n)), flipud(newL(top:m, left))'];
            label = searchArea(find(searchArea ~= -1, 1));
            if (~isempty(label))
                newL(noConnectedInd) = label;
                break
            end
            searchsize = searchsize + 1;
        end
    end
end

