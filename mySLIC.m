function [L, num, centerFeatures] = mySLIC(I, k, m) % S = sqrt(N / k)
% A instance of SLIC algroithom, which is for super pixel segemantation
% Input:
% I -- feature layers of the image to be segemanted 
% k -- the number of super pixels you want 
% m -- the weight of the distance in x-y space
% Output:
% L -- the labeled image of I
% num -- the actual num of super pixels, since the sample distance may not
% meet the end of width or height
% centerFeatures - the features of each cluster center 

tic
[M, N, ~] = size(I);
S = round(sqrt(M * N / k));
[Ns, Ms] = meshgrid(S : S : N, S : S : M);
Ms = reshape(Ms, [], 1); Ns = reshape(Ns, [], 1);
seeds = [Ms, Ns];
num = size(seeds, 1);
[FX, FY] = gradient(double(sum(I, 3)));
Gamp = sqrt(FX.^2 + FY.^2);
for i = 1 : size(seeds, 1)
    local = Gamp(seeds(i,1)-1 : seeds(i,1)+1, seeds(i,2) - 1 : seeds(i,2) + 1);
    [~, index] = min(local);    
    seeds(i, 1) = index(1) - 2 + seeds(i, 1);
    seeds(i, 2) = index(2) - 2 + seeds(i, 2);    
end
% init label image L
L = -1 * ones(M, N);
DisI = Inf * ones(M, N);
% m, n, other features
centerFeatures = zeros(size(seeds, 1), size(I, 3) + 2);
for i = 1 : size(seeds, 1)
    centerFeatures(i, 1:2) = seeds(i,:);
    centerFeatures(i, 3:end) = reshape(I(seeds(i,1), seeds(i, 2), :), 1, []);
end
times = 500;
for j = 1 : times
    isStable = true;
    for i = 1 : size(centerFeatures, 1)
        % search in neighboor
        mSearchRange = max(1, centerFeatures(i,1) - S) : min(centerFeatures(i,1) + S, M);
        nSearchRange = max(1, centerFeatures(i,2) - S) : min(centerFeatures(i,2) + S, N);
        localM = length(mSearchRange); localN = length(nSearchRange);
        localI = I(mSearchRange, nSearchRange, :);
        localDis = DisI(mSearchRange, nSearchRange);
        centerM = centerFeatures(i, 1); centerN = centerFeatures(i, 2);
        otherFeature = reshape(centerFeatures(i,3:end), 1, 1, []);
        % distance in otherfeature space
        newDislocal = sum((localI - repmat(otherFeature, localM, localN)).^2 , 3);       
        % distance in real space
        [Ncoor, Mcoor] = meshgrid(nSearchRange, mSearchRange);
        newDislocal = newDislocal + ((Mcoor - centerM).^2 + (Ncoor - centerN).^2) / S.^2 * m^2;
        newDislocal = sqrt(newDislocal);
        mask = newDislocal <  localDis;
        if (~isempty(find(mask)))
            isStable = false;
        end
        newL = L(mSearchRange, nSearchRange);
        newL(mask) = i;
        L(mSearchRange, nSearchRange) = newL;
        newDisBuffer = DisI(mSearchRange, nSearchRange);
        newDisBuffer(mask) = newDislocal(mask);
        DisI(mSearchRange, nSearchRange) = newDisBuffer;   
    end
    
    % update center
    for i = 1 : size(centerFeatures, 1)
        mSearchRange = max(1, centerFeatures(i,1) - S) : min(centerFeatures(i,1) + S, M);
        nSearchRange = max(1, centerFeatures(i,2) - S) : min(centerFeatures(i,2) + S, N);
        localM = length(mSearchRange); localN = length(nSearchRange);
        localI = I(mSearchRange, nSearchRange,:);
        localL = L(mSearchRange, nSearchRange,:);
        [Ncoor, Mcoor] = meshgrid(nSearchRange, mSearchRange);
        localFeature = zeros(localM * localN, size(I, 3) + 2);
        localFeature(:,1:2) = [reshape(Mcoor, [],1), reshape(Ncoor, [],1)];
        for channel = 1 : size(I, 3)
            localFeature(:, channel + 2) = reshape(localI(:,:,channel), [], 1);
        end
        % get means in cluster
        localFeature(find(localL ~= i),:) = [];
        meanFeature = mean(localFeature, 1);
        meanM = meanFeature(1); meanN = meanFeature(2);
        meanOther = reshape(meanFeature(3:end), 1, 1, []);
        % calulate the distance to mean
        % distance in otherfeature space
        toMeanDis = sum((localI - repmat(meanOther, localM, localN)).^2 , 3);     
        % distance in real space
        toMeanDis = toMeanDis + ((Mcoor - meanM).^2 + (Ncoor - meanN).^2) / S.^2 * m^2;
        toMeanDis = sqrt(toMeanDis);
        [~, pos] = min(toMeanDis(:));
        
        % update
        [m, n] = ind2sub(size(toMeanDis), pos);
        m = m(1); n = n(1);
        centerFeatures(i,:) = [Mcoor(m,n), Ncoor(m,n), reshape(localI(m,n,:),1,[])];
    end
    
    if (mod(j, 10) == 0)
        progress = j
    end
    if (isStable)
        break
    end
end
toc
end