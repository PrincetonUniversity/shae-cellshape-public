%% RELATED TO DVPSF
% select a folder with a bunch of psf stacks and separate into individual
% peaks and average them together

%% (1) select a set of folders with images of bright beads and find peak positions
% location to store all of the mat files with processed information from
% the stack
storeFolder = '/Volumes/BrattonData/psf20141027/averagedData';
 pixSize = 0.080; % xy pixel size in um
 xyPSFsize = 20; % radius in pixels to extract for averaging PSF
 zPSFsize = 40; % radius in pixels to extract for averaging PSF
    
%% (1a) Select folder with PSF stacks
topFolder = uigetdir();

dirList = dir(topFolder);
isDirIDX = [dirList.isdir];
isDirIDX(1:2) = 0;
dirList = dirList(isDirIDX);
    newsub=dir([topFolder filesep dirList(1).name]);

while isempty(strfind([newsub.name],'tif'))&&~isempty([newsub.name])
    subdir=dir([topFolder filesep dirList(1).name]);
    subdir(1:2)=[];
    if strfind(subdir(1).name,'DS')
        subdir(1)=[];
    end
    for ii=1:length(subdir)
        subdir(ii).name=[dirList(1).name filesep subdir(ii).name];
    end
    dirList=cat(1,dirList,subdir);
    dirList(1)=[];
    newsub=dir([topFolder filesep dirList(1).name]);
    
end

    




%% (1b) Find peaks, cut them out, move .mat files to communal location
% progressbar(0);
iStart = 1;
% iEnd = 1;
iEnd = length(dirList);
% iEnd = 2;
for iFolder = iStart:iEnd;
    % for iFolder = 1:32;
    % for iFolder = 4:4;
        progressbar((iFolder-iStart)/(iEnd-iStart));
    
    
    %     % smoothing/zeroing
    param.pedestal = 100; % pedestal value
    param.bpassNoise = 1;% noise lengthscale for 2D bandpass, pixels
    param.bpassFeature = 11; % feature lengthscale for 2D bandpass, pixels
    param.pixSize = pixSize; % xy pixel size in um
    param.kernelZ = 0.1; % Gaussian smooth along focal dimension in um
    
    %     % peak determination
    param.pkThresh = 20; % thresh for pkfnd
    param.pkSize = 5; % size for pkfnd
    
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
    param.dispPause = [];
    
    param.xyPSFsize = xyPSFsize; % radius in pixels to extract for averaging PSF
    param.zPSFsize = zPSFsize; % radius in pixels to extract for averaging PSF
    
    param.saveFile = true;
    param.isInterp3D = false;
    
    % for iFolder = 1;
    find3DPeaksAtSurface3(fullfile(topFolder,dirList(iFolder).name),param);
    matFile = dir([fullfile(topFolder,dirList(iFolder).name),filesep,'*.mat']);
    movefile(fullfile(topFolder,dirList(iFolder).name,matFile(1).name),...
        storeFolder);
end
progressbar(1);
%% (2) given a set of .mat files with peak positions, process to average PSF
%% (2a) select files

[FileName,PathName] = uigetfile('','MultiSelect','on');
if not(iscell(FileName))
    FileName = {FileName};
end



%% (2b) group PSFs into sections, get set up
zPSFsectionCenter = -5:0.5:1;
param.xyPSFsize = xyPSFsize; % radius in pixels to extract for averaging PSF
param.zPSFsize = zPSFsize; % radius in pixels to extract for averaging PSF

% set up array
averagedPSF = struct();
for kkSlice = 1:size(zPSFsectionCenter,2)
    averagedPSF(kkSlice).addedPSF = zeros(param.xyPSFsize*2+1,param.xyPSFsize*2+1,param.zPSFsize*2+1);
    averagedPSF(kkSlice).sumMask=zeros(param.xyPSFsize*2+1,param.xyPSFsize*2+1,param.zPSFsize*2+1);
    averagedPSF(kkSlice).pkIntensityArray = [];
    averagedPSF(kkSlice).pkBkgdArray = [];
    averagedPSF(kkSlice).zStep = [];
end
% plane
planeFunc = @(a,x) a(1)+a(2).*x(:,1)+a(3).*x(:,2);

%% (2c) load files individually and estimate surface
for iiFile = 1:size(FileName,2)
    % for iiFile = 1;
    % load file
    load(fullfile(PathName,FileName{iiFile}));
    % make sure there are a reasonable number
    if size(psfExample,2)<5
 
        continue
    end
    %% (2c1) Fit all points to a single plane
    xyzSurf = nan(size(psfExample,2),3);
    for jjPeak = 1:size(psfExample,2)
        xyzSurf(jjPeak,1) = psfExample(jjPeak).xCentPix;
        xyzSurf(jjPeak,2) = psfExample(jjPeak).yCentPix;
        xyzSurf(jjPeak,3) = psfExample(jjPeak).bestEst;
    end
    
    % fit to plane
    [bestPlaneOut,resnorm,residual] = lsqcurvefit(planeFunc,[1,1,1],...
        xyzSurf(:,1:2),xyzSurf(:,3));
    
    %% (2c2) group/fit residuals to two gaussians
    try
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
    catch
    end
    
    %% (2c3) group/fit residuals to two gaussians, again
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
    
    %% (2c4) convert z positions into relative to surface positions
    xyzSurf = nan(size(psfExample,2),3);
    for jjPeak = 1:size(psfExample,2)
        xyzSurf(jjPeak,1) = psfExample(jjPeak).xCentPix;
        xyzSurf(jjPeak,2) = psfExample(jjPeak).yCentPix;
        xyzSurf(jjPeak,3) = psfExample(jjPeak).bestEst;
    end
    
    xyzSurf(:,3) = xyzSurf(:,3)-planeFunc(bestPlaneOut,xyzSurf(:,1:2));
    
    
    
    %% (2c5) determine which section each peak falls within
    [n,zPeakSectionID] = histc(xyzSurf(:,3),zPSFsectionCenter);
    
    %% (2c6) recast image as doubles, normalize and add to others
    for jjPeak = 1:size(xyzSurf,1)
        tempPSF = psfExample(jjPeak).psfCutOut;
        tempPSF = double(tempPSF);
        tempPSF = tempPSF./nanmax(tempPSF(:));
        tempMask=~isnan(tempPSF);
        tempPSF(isnan(tempPSF))=0;
        
        
        % append data to appropriate z section/slice
        if zPeakSectionID(jjPeak)>0 % only use it if it is in the appropriate range
            if isempty(averagedPSF(zPeakSectionID(jjPeak)).zStep)
                if isfield(psfExample(jjPeak),'zStep')
                    
                averagedPSF(zPeakSectionID(jjPeak)).zStep = psfExample(jjPeak).zStep;
                end
            else
                if isfield(psfExample(jjPeak),'zStep')
                    if averagedPSF(zPeakSectionID(jjPeak)).zStep ~= psfExample(jjPeak).zStep
                        warning('processPSFStacks:differentZStep','Different step size along focal dimension');
                        continue
                    end
                end
            end
            
            averagedPSF(zPeakSectionID(jjPeak)).addedPSF = averagedPSF(zPeakSectionID(jjPeak)).addedPSF+tempPSF;

            averagedPSF(zPeakSectionID(jjPeak)).pkIntensityArray = cat(2,...
                averagedPSF(zPeakSectionID(jjPeak)).pkIntensityArray,...
                psfExample(jjPeak).psfPeak);
            averagedPSF(zPeakSectionID(jjPeak)).sumMask=averagedPSF(zPeakSectionID(jjPeak)).sumMask +tempMask;
            
            averagedPSF(zPeakSectionID(jjPeak)).pkBkgdArray = cat(2,...
                averagedPSF(zPeakSectionID(jjPeak)).pkBkgdArray,...
                psfExample(jjPeak).psfBkgd);
        end % appending peak data to proper z section/slice
    end % loop over peaks
end % loop over .mat files

%% (3)
% remember to save the mat file at some point

%% (4) query total number of spots averaged together in each slice
nDotsAvgKKSlice = [];
medianInte = [];
meanInte = [];
quantileA = [];
quantileB = [];
distIntenZ = [];

for kkSlice = 1:length(zPSFsectionCenter)-1
    nDotsAvgKKSlice(kkSlice) = length(averagedPSF(kkSlice).pkIntensityArray);
    meanInte(kkSlice) = expfit(averagedPSF(kkSlice).pkIntensityArray(:));
end

%% um, don't remember this section
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
%% plot the mean intensity for brightness per slice
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

% %% load file and estimate surface
% addedPSF = zeros(size(psfExample(1).psfCutOut));
% counter = 0;
% % for iiFile = 1:size(FileName,2)
%     for iiFile = 1:2s;
%     % load file
%     load(fullfile(PathName,FileName{iiFile}));
%     
%     for jjPeak = 1:size(psfExample,2)
%         tempPSF = psfExample(jjPeak).psfCutOut;
%         tempPSF = double(tempPSF);
%         tempPSF = tempPSF./max(tempPSF(:));
%         addedPSF = addedPSF+tempPSF;
%         counter = counter+1;
%     end % loop over peaks
% end % loop over .mat files

