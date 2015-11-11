%% RELATED TO DVPSF
% select a folder with a bunch of psf stacks and separate into individual
% peaks and average them together

%% Select folder with PSF stacks
topFolder = uigetdir();

dirList = dir(topFolder);
isDirIDX = [dirList.isdir];
isDirIDX(1:2) = 0;
dirList = dirList(isDirIDX);
%% Find peaks, cut them out, move .mat files to communal location
% progressbar(0);
iStart = 1;
% iEnd = 1;
iEnd = length(dirList);
iEnd = 4;
for iFolder = iStart:iEnd;
    % for iFolder = 1:32;
    % for iFolder = 4:4;
    %     progressbar((iFolder-iStart)/(iEnd-iStart));
    
    
    %     % smoothing/zeroing
    param.pedestal = 100; % pedestal value
    param.bpassNoise = 1;% noise lengthscale for 2D bandpass, pixels
    param.bpassFeature = 11; % feature lengthscale for 2D bandpass, pixels
    param.pixSize = 0.085; % xy pixel size in um
    param.kernelZ = 0.1; % Gaussian smooth along focal dimension in um
    
    %     % peak determination
    param.pkThresh = 3; % thresh for pkfnd
    param.pkSize = 4; % size for pkfnd
    
    %     % finding z position
    param.zRadius = 1; % radius, in microns to fit intensity(z) to polynomial
    param.nPoly = 4; % order of polynomial to fit intensity(z)
    
    % radius, in pixels, of box to draw around centroid of thresholded
    % image to find peak pixel
    param.radiusBox = 3;
    % radius, in pixels, of box to integrate around peak pixel for determining intensity(z)
    param.radiusIntegral = 7;
    
    %     % properties of intensity(z)
    param.zIntensityMin = 100; % the integral of the peak must be at least this in one plane
    param.tempA = 2; % total number of local peaks (countsA+countsB);
    % % false/noise peaks
    param.tempB1 = 0.35;
    param.tempB2 = 0.6;
    param.tempA1 = 0.15;
    param.tempA2 = 0.25;
    % % xy center and size determination
    param.avgZForXY = 1; % radius along z to average for determination of XY centroid
    
    param.similarBestEst = 0.750; % um different between z in focus estimates allowed
    
    % % 'intensity' mode, now in options
    %     param.intensityMode = 'Brenner'; % use the Brenner gradient as a way of determining in focus
    %           param.intensityMode = 'Sum'; % use the integral in a region around
    %       the center in xy as a metric for in focus
    
    % % various display modes
    param.displayIntensity = 0; % set to true to display intensity graphs
    param.dispPause = 0;
    
    param.xyPSFsize = 25; % radius in pixels to extract for averaging PSF
    param.zPSFsize = 25; % radius in pixels to extract for averaging PSF
    
    
    % for iFolder = 1;
    [savedPeaks] = find3DPeaksAtSurface(fullfile(topFolder,dirList(iFolder).name),param);
    matFile = dir([fullfile(topFolder,dirList(iFolder).name),filesep,'*.mat']);
    movefile(fullfile(topeFolder,dirList(iFolder).name,matFile(1).name),...
        ['/Users/benbratton/Documents/MATLAB/misc matlab data/slantA/']);




end
progressbar(1);

%% given a .mat files with peak positions, process to average PSF

[FileName,PathName] = uigetfile('','MultiSelect','on');
if not(iscell(FileName))
    FileName = {FileName};
end



%% group PSFs by 100 nm sections
zPSFsectionCenter = -5:0.05:1;
param.xyPSFsize = 25; % radius in pixels to extract for averaging PSF
param.zPSFsize = 15; % radius in pixels to extract for averaging PSF

% set up array
averagedPSF = struct();
for kkSlice = 1:size(zPSFsectionCenter,2)
    averagedPSF(kkSlice).addedPSF = zeros(param.xyPSFsize*2+1,param.xyPSFsize*2+1,param.zPSFsize*2+1);
    averagedPSF(kkSlice).pkIntensityArray = [];
    averagedPSF(kkSlice).pkBkgdArray = [];
end
%% plane
planeFunc = @(a,x) a(1)+a(2).*x(:,1)+a(3).*x(:,2);

%% load file and estimate surface
for iiFile = 1:size(FileName,2)
    % for iiFile = 1;
    % load file
    load(fullfile(PathName,FileName{iiFile}));
    % make sure there are a reasonable number
    if size(psfExample,2)<20
        continue
    end
    %%
    xyzSurf = nan(size(psfExample,2),3);
    for jjPeak = 1:size(psfExample,2)
        xyzSurf(jjPeak,1) = psfExample(jjPeak).xCentPix;
        xyzSurf(jjPeak,2) = psfExample(jjPeak).yCentPix;
        xyzSurf(jjPeak,3) = psfExample(jjPeak).bestEst;
    end
    
    % fit to plane
    [bestPlaneOut,resnorm,residual] = lsqcurvefit(planeFunc,[1,1,1],...
        xyzSurf(:,1:2),xyzSurf(:,3));
    
    % fit residuals to two gaussians
    gmdist1 = gmdistribution.fit(residual,2,'Replicates',10);
    % which component has small variance?
    [~,idxComp] = min(gmdist1.Sigma,[],3);
    % which peaks come primarily from that distribution?
    residPdf1 = pdf('norm',residual,gmdist1.mu(idxComp),gmdist1.Sigma(1,1,idxComp));
    if idxComp == 1
        residPdf2 = pdf('norm',residual,gmdist1.mu(2),gmdist1.Sigma(1,1,2));
    else
        residPdf2 = pdf('norm',residual,gmdist1.mu(1),gmdist1.Sigma(1,1,1));
    end
    
    idxKeepPlane = residPdf1>residPdf2;
    xyzSurf = xyzSurf(idxKeepPlane,:);
    
    [bestPlaneOut,resnorm,residual] = lsqcurvefit(planeFunc,[1,1,1],...
        xyzSurf(:,1:2),xyzSurf(:,3));
    % fit residuals to two gaussians
    try
        gmdist1 = gmdistribution.fit(residual,2,'Replicats',10);
        % which component has big intensity or lots of peaks?
        [~,idxComp] = max(gmdist1.PComponents,[],2);
        % which peaks come primarily from that distribution?
        residPdf1 = pdf('norm',residual,gmdist1.mu(idxComp),gmdist1.Sigma(1,1,idxComp));
        if idxComp == 1
            residPdf2 = pdf('norm',residual,gmdist1.mu(2),gmdist1.Sigma(1,1,2));
        else
            residPdf2 = pdf('norm',residual,gmdist1.mu(1),gmdist1.Sigma(1,1,1));
        end
        
        idxKeepPlane = residPdf1>residPdf2;
        xyzSurf = xyzSurf(idxKeepPlane,:);
    catch ME % keep everyone left
    end
    
    
    [bestPlaneOut,resnorm,residual] = lsqcurvefit(planeFunc,[1,1,1],...
        xyzSurf(:,1:2),xyzSurf(:,3));
    
    % convert z positions into relative to surface positions
    xyzSurf = nan(size(psfExample,2),3);
    for jjPeak = 1:size(psfExample,2)
        xyzSurf(jjPeak,1) = psfExample(jjPeak).xCentPix;
        xyzSurf(jjPeak,2) = psfExample(jjPeak).yCentPix;
        xyzSurf(jjPeak,3) = psfExample(jjPeak).bestEst;
    end
    
    xyzSurf(:,3) = xyzSurf(:,3)-planeFunc(bestPlaneOut,xyzSurf(:,1:2));
    
    
    
    %%
    % determine which section each peak falls within
    [n,zPeakSectionID] = histc(xyzSurf(:,3),zPSFsectionCenter);
    
    %% recast as doubles and add to others
    for jjPeak = 1:size(xyzSurf,1)
        tempPSF = psfExample(jjPeak).psfCutOut;
        tempPSF = double(tempPSF);
        tempPSF = tempPSF./max(tempPSF(:));
        % append data to appropriate z section/slice
        if zPeakSectionID(jjPeak)>0 % only use it if it is in the appropriate range
            averagedPSF(zPeakSectionID(jjPeak)).addedPSF = averagedPSF(zPeakSectionID(jjPeak)).addedPSF+tempPSF;
            averagedPSF(zPeakSectionID(jjPeak)).pkIntensityArray = cat(2,...
                averagedPSF(zPeakSectionID(jjPeak)).pkIntensityArray,...
                psfExample(jjPeak).psfPeak);
            averagedPSF(zPeakSectionID(jjPeak)).pkBkgdArray = cat(2,...
                averagedPSF(zPeakSectionID(jjPeak)).pkBkgdArray,...
                psfExample(jjPeak).psfBkgd);
        end % appending peak data to proper z section/slice
    end % loop over peaks
end % loop over .mat files
%%

%% query total number of qDots averaged together in each slice
nDotsAvgKKSlice = [];
medianInte = [];
meanInte = [];
quantileA = [];
quantileB = [];
% binsInten = 0:45:750;
distIntenZ = [];

for kkSlice = 1:length(zPSFsectionCenter)-1
    nDotsAvgKKSlice(kkSlice) = length(averagedPSF(kkSlice).pkIntensityArray);
    %    medianInte(kkSlice) = median(averagedPSF(kkSlice).pkIntensityArray(:));
    meanInte(kkSlice) = expfit(averagedPSF(kkSlice).pkIntensityArray(:));
    %    quantileA(kkSlice) = quantile(averagedPSF(kkSlice).pkIntensityArray(:),0.75);
    %    quantileA(kkSlice) = quantile(averagedPSF(kkSlice).pkBkgdArray(:),0.75);
    %    tempCounts = histc(averagedPSF(kkSlice).pkIntensityArray(:),binsInten);
    %    distIntenZ = cat(2,distIntenZ,tempCounts./sum(tempCounts));
end
%
% imHandle = imshow(distIntenZ,[],'initialMagnification','fit');
% set(get(imHandle,'Parent'),'XTickLabel',zPSFsectionCenter(get(get(imHandle,'Parent'),'XTick')));
% set(get(imHandle,'Parent'),'YTickLabel',binsInten(get(get(imHandle,'Parent'),'YTick')));
% set(get(imHandle,'Parent'),'YDir','normal');
% set(get(imHandle,'Parent'),'Visible','On');
% ylabel('peak intensity');
% xlabel('z position');
% cbarHandle = colorbar();
% xlabel(cbarHandle,'p(Int|z)');
% % scatter(zPSFsectionCenter(1:end-1),quantileA);
% % ylim([0,750]);
avgInt = [];
for ii = 1:120/5
    intensityArray = [];
    for jj = 1:5
        intensityArray = cat(1,...
            averagedPSF(5*(ii-1)+jj).pkIntensityArray(:));
    end
    avgInt(ii,1) = zPSFsectionCenter(5*(ii-1)+3);
    [avgInt(ii,2),avgInt(ii,3:4)] = expfit(intensityArray);
end
%%
clf;
hold on;
scatter(zPSFsectionCenter(1:120),meanInte,'kx');
% plot(sgolayfilt(meanInte,3,15),'r');
plot(zPSFsectionCenter(1:120), sgolayfilt(meanInte,3,45),'b');
% scatter(avgInt(:,1),avgInt(:,2),'ro');
% scatter(avgInt(:,1),avgInt(:,3),'rs');
% scatter(avgInt(:,1),avgInt(:,4),'rs');
ylim([0,1000]);
xlabel('height in sample \mum');
ylabel('average peak pixel intensity');
% scatter(zPSFsectionCenter(1:120),nDotsAvgKKSlice);

%% %% %%
%% load file and estimate surface
addedPSF = zeros(size(psfExample(1).psfCutOut));
counter = 0;
% for iiFile = 1:size(FileName,2)
    for iiFile = 1:2s;
    % load file
    load(fullfile(PathName,FileName{iiFile}));
    
    for jjPeak = 1:size(psfExample,2)
        tempPSF = psfExample(jjPeak).psfCutOut;
        tempPSF = double(tempPSF);
        tempPSF = tempPSF./max(tempPSF(:));
        addedPSF = addedPSF+tempPSF;
        counter = counter+1;
    end % loop over peaks
end % loop over .mat files


%% group multiple files together 

% select a few files
whichFiles = uigetfile('*.mat','MultiSelect','on');
if not(iscell(whichFiles))
    whichFiles = {whichFiles};
end
savedPeaks2 = [];
for iFile = 1:size(whichFiles,2);
   load(whichFiles{iFile},'savedPeaks','firstStagePos','emFilter');
   switch emFilter
       case '534/30 nm (2)'
           savedPeaks(:,18) = 2;
       case '641/75 nm (3)'
           savedPeaks(:,18) = 3;
       otherwise
           savedPeaks(:,18) = nan;
   end
   savedPeaks(:,1) = savedPeaks(:,1)/1000+firstStagePos(1);
   savedPeaks(:,2) = savedPeaks(:,2)/1000+firstStagePos(2);
   savedPeaks(:,17) = savedPeaks(:,17)*(10^iFile);
   savedPeaks2 = cat(1,savedPeaks2,savedPeaks);
end

%%
hist(savedPeaks(:,3),100);
savedPeaks2(savedPeaks2(:,3)==0,:)=[];
newVar = [savedPeaks2(:,3),log(savedPeaks2(:,17))/6];
% newVar(newVar(:,2)<1.5,:) = [];

%% remove asymetric peaks
keepIdx = (savedPeaks2(:,14)./savedPeaks2(:,15))<1.15;
savedPeaks3 = savedPeaks2(keepIdx,:);
newVar = [savedPeaks3(:,3),log(savedPeaks3(:,17))/6];

[IDX,C] = kmeans(newVar,7,'replicates',50);

clf
hold on;
colrs = jet(max(IDX));
for ii = 1:max(IDX);
% for ii = 3;
scatter(newVar(IDX==ii,1),(newVar(IDX==ii,2)),15,colrs(ii,:),'o');
scatter(C(ii,1),C(ii,2),35,colrs(ii,:),'x','LineWidth',2);
end
axis equal;
figure(gcf);

%%
figure(gcf);
clf;
idx = savedPeaks2(:,18);
scatter3(savedPeaks2(idx==2,1),savedPeaks2(idx==2,2),savedPeaks2(idx==2,3),'gx');
hold on;
scatter3(savedPeaks2(idx==3,1),savedPeaks2(idx==3,2),savedPeaks2(idx==3,3),'rx');

%%

planeFunc = @(a,x) a(1)+a(2).*x(:,1)+a(3).*x(:,2);
    % fit to plane
    [bestPlaneOut,resnorm,residual] = lsqcurvefit(planeFunc,[1,1,1],...
        savedPeaks2(IDX==ii,1:2),savedPeaks2(IDX==ii,3));

    new3 = savedPeaks2;
    new3(:,3) = new3(:,3)-planeFunc(bestPlaneOut,new3(:,[1,2]));
    
    scatter(new3(:,3),log(new3(:,17)));