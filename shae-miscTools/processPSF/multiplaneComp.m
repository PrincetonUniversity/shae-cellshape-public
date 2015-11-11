function residuals = multiplaneComp(a,x,zIDList,zPrimeList,dvpsfStruct);
% this function is just the residuals, the minimization will square and
% sum

% check to make sure that there are same number of zPrime planes in
% dvpsfStruct and the list of their coordinates
%% utilize the parameters in a to calculate single planes

if size(zPrimeList,1)~=size(dvpsfStruct,1)
    error('multiplaneComp:sizeMismatch:zPrime','improper number of depth varying sections');
end

if max(zIDList)>size(dvpsfStruct(1).addedPSF,3)
    error('multiplaneComp:sizeMismatch:zID','improper focal section selection');
end

if min(zIDList)<1
    error('multiplaneComp:sizeMismatch:zID','improper focal section selection');
end

%% loop over places in sample
residuals = [];
for iiZPrime = 1:size(zPrimeList,2)
    % extract psf
    tempPSF = dvpsfStruct(iiZPrime).addedPSF;
% assume peak intensity is unaffected, reasonable for first 1-2 um into
% sample
    %     meanIntensity = mean(dvpsfStruct(iiZPrime).pkIntensityArray);
    
    % just look at the planes of interest
    tempPSF = tempPSF(:,:,zIDList);
    
    % renormalize to a peak voxel of 1
    tempPSF = tempPSF./max(tempPSF(:));
    
    % calculate multiplane
    [out,bSave] = multiPlane(a,...
        x,...
        zIDList-mean(zIDList),...
        zPrimeList(iiZPrime),...
        [size(tempPSF,1),size(tempPSF,2)]);
    
    assignin('base','bSave',bSave);
    % calculate residuals
    
    residuals(:,iiZPrime) = tempPSF(:) - out(:);
    
end % iiZprime, place in sample

% return residuals

