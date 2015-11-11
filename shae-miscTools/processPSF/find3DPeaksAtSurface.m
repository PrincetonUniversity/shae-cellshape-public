function [savedPeaks] = find3DPeaksAtSurface(folderName,param)
if nargin<1
    %% open an image
    folderName = ['/Users/benbratton/Documents/data from theBlock/20120726_colorShiftIssues/20120726_163914'];
    % list files to read
end

%% set up peak finding parameters
if nargin<1
    
    %     % smoothing/zeroing
    param.pedestal = 100; % pedestal value
    param.bpassNoise = 1;% noise lengthscale for 2D bandpass, pixels
    param.bpassFeature = 11; % feature lengthscale for 2D bandpass, pixels
    param.pixSize = 0.075; % xy pixel size in um
    param.kernelZ = 0.1; % Gaussian smooth along focal dimension in um
    
    %     % peak determination
    param.pkThresh = 7; % thresh for pkfnd
    param.pkSize = 8; % size for pkfnd
    
    %     % finding z position
    param.zRadius = 1; % radius, in microns to fit intensity(z) to polynomial
    param.nPoly = 4; % order of polynomial to fit intensity(z)
    
    % radius, in pixels, of box to draw around centroid of thresholded
    % image to find peak pixel
    param.radiusBox = 3;
    % radius, in pixels, of box to integrate around peak pixel for determining intensity(z)
    param.radiusIntegral = 3;
    
    %     % properties of intensity(z)
    param.zIntensityMin = 200; % the integral of the peak must be at least this in one plane
    param.tempA = 2; % total number of local peaks (countsA+countsB);
    % % false/noise peaks
    param.tempB1 = 0.35;
    param.tempB2 = 0.6;
    param.tempA1 = 0.15;
    param.tempA2 = 0.25;
    % % xy center and size determination
    param.avgZForXY = 1; % radius along z to average for determination of XY centroid
    
    param.similarBestEst = 0.250; % um different between z in focus estimates allowed
    
    % % 'intensity' mode, now in options
    %     param.intensityMode = 'Brenner'; % use the Brenner gradient as a way of determining in focus
    %           param.intensityMode = 'Sum'; % use the integral in a region around
    %       the center in xy as a metric for in focus
    
    % % various display modes
    param.displayIntensity = 0; % set to true to display intensity graphs
    param.dispPause = 0;
    
    param.xyPSFsize = 25; % radius in pixels to extract for averaging PSF
    param.zPSFsize = 15; % radius in pixels to extract for averaging PSF

    param.saveFile = true;
end

%% prepare to read in stack
fileList = dir(fullfile(folderName,'*.tif*'));

% read in each file
imageInfo = {[]};
tempImage = [] ;
zPos = [];

%% read in stack
for iFile = 1:length(fileList);
    % note: this construction can return the original array without any
    % filtering
    tempImage(:,:,iFile) = bpass(double(imread(...
        fullfile(folderName,fileList(iFile).name))),...
        0,0);
    %     pkFindImage(:,:,iFile) = bpass(double(imread(...
    %         fullfile(folderName,fileList(iFile).name))),...
    %         1,5);
    
    
    imageInfo{iFile} = parseTags(  ...
        fullfile(folderName,fileList(iFile).name));
    if not(isempty(fieldnames(imageInfo{iFile})))
        zPos(iFile) = imageInfo{iFile}.posZ;
    end
end

if isempty(zPos)
    %sliceSize = inputdlg('Spacing between planes in um?');
    sliceSize='0.1';
    zPos = (1:size(tempImage,3)).*str2double(sliceSize);
end

% the step size between planes in stage microns
zStep = zPos(2)-zPos(1);

%% remove the pedestal
tempImage = tempImage-param.pedestal;

%% smooth each plane in xy to smooth out pixel to pixel variations and uneven background
tempImage2 = nan(size(tempImage));
for iPlane = 1:size(tempImage,3)
    tempImage2(:,:,iPlane) = bpass(tempImage(:,:,iPlane),param.bpassNoise,param.bpassFeature);
    
end


% smooth a second time in z
z2 = -3:3;
kernel = exp(-(z2.^2)/(2*(param.kernelZ/(zPos(2)-zPos(1)))^2));
kernel = kernel./sum(kernel(:));
kernel = reshape(kernel,[1,1,7]);

tempImage2 = imfilter(tempImage2,reshape(kernel,[1,1,7]));

%% pkfind on std image
peaksOut = pkfnd(std(tempImage2,0,3),param.pkThresh,param.pkSize);

%% go through possible peaks, centered on peak pixel and calculate intensity as a function of Z
peaks2 = nan(size(peaksOut,1),11);
psfExample = [];
% for iPeak = 261;
psfCounter = 0;

for iPeak = 1:size(peaksOut,1);
    if peaksOut(iPeak,1)<param.radiusIntegral +1 || peaksOut(iPeak,1)>size(tempImage,1)-param.radiusIntegral  ...
            || peaksOut(iPeak,2)<param.radiusIntegral +1 || peaksOut(iPeak,2)>size(tempImage,2)-param.radiusIntegral ;
        continue
    end
    % xPeak, yPeak of max pixel;
    xPeak = peaksOut(iPeak,1)*param.pixSize*1000;
    yPeak = peaksOut(iPeak,2)*param.pixSize*1000;
    peaks2(iPeak,1:2) = [xPeak,yPeak];
    
    % cut out near center of peak for calculating zPosition
    pkCutOut = tempImage2(peaksOut(iPeak,2)-param.radiusIntegral:peaksOut(iPeak,2)+param.radiusIntegral ,...
        peaksOut(iPeak,1)-param.radiusIntegral :peaksOut(iPeak,1)+param.radiusIntegral ,...
        :);
    
    % run peak estimator with sum intensity
    options.intensityMode = 'Sum';
    [bestEstSum,outputSum] = ZPeakEstimator(pkCutOut,param,zPos,options);
    
    % run peak estimator with Brenner gradient
    options.intensityMode = 'Brenner';
    [bestEstBren,outputBren] = ZPeakEstimator(pkCutOut,param,zPos,options);
    
    % store results
    peaks2(iPeak,5) = bestEstBren;
    peaks2(iPeak,6) = outputBren.isKept;
    peaks2(iPeak,7) = outputBren.fwhmZ;
    peaks2(iPeak,8) = bestEstSum;
    peaks2(iPeak,9) = outputSum.isKept;
    peaks2(iPeak,10) = outputSum.fwhmZ;
    
    % if both peaks are accepted and they give similar estimates of z
    % position, keep them
    if and(and(abs(bestEstBren-bestEstSum)<param.similarBestEst,...
            outputSum.isKept),outputBren.isKept);
        
        % store results
        peaks2(iPeak,4) = 1;
        
        % extract a region around peak in xy and pass to subroutine
        bestEst = ((bestEstBren+bestEstSum)/2);
        bestEstP = round((bestEst-zPos(1))/(zPos(2)-zPos(1))+1);
        
        
        subImageIn =     tempImage2(peaksOut(iPeak,2)-param.pkSize:peaksOut(iPeak,2)+param.pkSize,...
            peaksOut(iPeak,1)-param.pkSize:peaksOut(iPeak,1)+param.pkSize,...
            bestEstP-param.avgZForXY:bestEstP+param.avgZForXY);
        subImageIn = mean(subImageIn,3);
        try
        [xCent,yCent,majAxis,minAxis,orientation] = ...
            xyCentSize(subImageIn);
        catch ME
            xCent = nan;
            yCent = nan;
            majAxis = nan;
            minAxis = nan;
            orientation = nan;
        end
            
        %         interactive = true;
        %         aab = gaussFit(subImageIn,[param.pkSize,param.pkSize],param.pkSize*2+1,interactive)
        
        % store centroid results
        xCntrd = 1000*param.pixSize*( xCent+peaksOut(iPeak,1)-param.pkSize-1);
        yCntrd = 1000*param.pixSize*( yCent+peaksOut(iPeak,2)-param.pkSize-1);
        
        % calculate total intensity/brightness
        subImageIn = tempImage2(peaksOut(iPeak,2)-param.pkSize:peaksOut(iPeak,2)+param.pkSize,...
            peaksOut(iPeak,1)-param.pkSize:peaksOut(iPeak,1)+param.pkSize,...
            min(outputSum.fwhmZ1,outputSum.fwhmZ2):max(outputSum.fwhmZ1,outputSum.fwhmZ2));
        
        
        [maxVal,maxId] = max(subImageIn(:));
        [~,~,tempZ] = ind2sub(size(subImageIn),maxId);
        maskXY = subImageIn(:,:,tempZ)>=maxVal.*0.5;
        totalIntensity = sum(sum(sum(bsxfun(@times,subImageIn,single(maskXY)),1),2),3);
        
        % store results
        peaks2(iPeak,3) = bestEst;
        %         peaks2(iPeak,11) = ?
        peaks2(iPeak,12) = xCntrd;
        peaks2(iPeak,13) = yCntrd;
        peaks2(iPeak,14) = majAxis;
        peaks2(iPeak,15) = minAxis;
        peaks2(iPeak,16) = orientation;
        peaks2(iPeak,17) = totalIntensity;
        
        % % calculate relative psf image
        
        % % make sure that the center is not near the edge
        %         psfExample will be the cell array to hold psf's
        newCentPixX = round(xCntrd/(1000*param.pixSize));
        newCentPixY = round(yCntrd/(1000*param.pixSize));
        
        
        % outside the range necessary for a full averaging
        if newCentPixX<param.xyPSFsize+1 || newCentPixX> size(tempImage,2)-param.xyPSFsize || ...
                newCentPixY<param.xyPSFsize+1 || newCentPixY> size(tempImage,1)-param.xyPSFsize || ...
                bestEstP<param.zPSFsize+1 || bestEstP> size(tempImage,3)-param.zPSFsize
        else
            
            
            % store some results
            psfCounter = psfCounter +1;
            psfExample(psfCounter).xCentPix = newCentPixX;
            psfExample(psfCounter).yCentPix = newCentPixY;
            psfExample(psfCounter).zCentPix = bestEstP;
            psfExample(psfCounter).bestEst  = bestEst;
            psfExample(psfCounter).pkID  = iPeak;
            
            try
            % extract image
            psfExample(psfCounter).psfCutOut = tempImage(newCentPixY-param.xyPSFsize:newCentPixY+param.xyPSFsize ,...
                newCentPixX-param.xyPSFsize:newCentPixX+param.xyPSFsize ,...
                bestEstP-param.zPSFsize:bestEstP+param.zPSFsize);
            
            % normalize to peak intensity
            psfExample(psfCounter).psfBkgd = median(psfExample(psfCounter).psfCutOut(:));
            psfExample(psfCounter).psfCutOut = psfExample(psfCounter).psfCutOut - psfExample(psfCounter).psfBkgd;
            psfExample(psfCounter).psfPeak = max(psfExample(psfCounter).psfCutOut(:));
            
            psfExample(psfCounter).psfCutOut = psfExample(psfCounter).psfCutOut./psfExample(psfCounter).psfPeak;
         psfExample(psfCounter).zStep = abs(zStep);
         
         if sign(zStep)<0
             psfExample(psfCounter).psfCutOut =flipdim(psfExample(psfCounter).psfCutOut,3);
         else
             
         end
            %             psfExample(psfCounter).psfCutOut = uint16(...
%                 psfExample(psfCounter).psfCutOut*((2^16)-1));
            psfExample(psfCounter).savedPeaksIndx = iPeak;
            catch ME
            end
        end
        
        
    else
        % store results
        peaks2(iPeak,4) = 0;
    end
    
    
    
    
    
end

%% remove peaks that aren't peaks
savedPeaks = peaks2;
% savedPeaks(savedPeaks(:,4)==0,:) = [];
savedPeaks = unique(savedPeaks,'rows');

%% store results
%%
saveFile = fullfile(folderName,[transcodeDate(now),'.mat']);


% emission filter is the filter used
% folderName is the path to the folder analyzed
% saved peaks is the xy/z coordinates of the peaks
%   % might be nice to have peak/integral intensities?
%   % might be nice to have fwhm in xy/z?
% param is the structure of various parameters used to find the peaks
% firstStagePos is the stage position of the first slice
% stagePosDif is direction and size of stage motion
try
    firstStagePos = [imageInfo{1}.posX,imageInfo{1}.posY,imageInfo{1}.posZ];
    stagePosDif = [imageInfo{2}.posX,imageInfo{2}.posY,imageInfo{2}.posZ]-firstStagePos;
catch
    firstStagePos = [0,0,0];
    stagePosDif = [0,0,sliceSize];
end

try
    emFilter = imageInfo{1}.filterID;
catch
    emFilter = '';
end

disp(folderName);   
if param.saveFile
save(saveFile,'folderName','savedPeaks','param','firstStagePos','stagePosDif','psfExample','emFilter');
disp(saveFile);
disp([num2str(psfCounter),' spots']);
elseif nargout<1
    warning('find3DPeaksAtSurface:paramMismatch:noOutput','Results will not be saved nor returned to the user');
end
% disp(thresh);
disp([]);
