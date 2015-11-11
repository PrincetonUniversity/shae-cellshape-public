%% (2) given a set of .mat files with peak positions, process to average PSF
%% (2a) select files

[FileName,PathName] = uigetfile('','MultiSelect','on');
if not(iscell(FileName))
    FileName = {FileName};
end


%% load in one at a time
% for iFile = 1:length(FileName)
surfacePSF = zeros(61,61,121);
surfacePSFcnts = zeros(61,61,121);
nBeads = 0;
for iFile = 1:36;
    progressbar((iFile-1)/135);
    beadStruct = load(fullfile(PathName, FileName{iFile}));
    listOfBeads = zeros(size(beadStruct.savedPeaks,1),4);
    for jBead = 1:size(beadStruct.psfExample,2)
       listOfBeads(beadStruct.psfExample(jBead).savedPeaksIndx,1) = jBead;
    end


beadPos = beadStruct.savedPeaks;
beadPos(not(listOfBeads(:,1)),:) = [];
listOfBeads(not(listOfBeads(:,1)),:) = [];
listOfBeads(:,2)  = beadPos(:,1)/1000;
listOfBeads(:,3)  = beadPos(:,2)/1000;
listOfBeads(:,4)  = beadPos(:,3);



rejectBeads = abs(zscore(listOfBeads(:,4)))>1.5;
%% find triangulation, shortest edges and remove
dt1 = delaunayTriangulation(listOfBeads(:,2),listOfBeads(:,3));

edgesList = dt1.edges;

edgeLengths = sqrt((listOfBeads(edgesList(:,1),1)-listOfBeads(edgesList(:,2),1)).^2+...
    (listOfBeads(edgesList(:,1),2)-listOfBeads(edgesList(:,2),2)).^2);
lengthCutoff = beadStruct.param.xyPSFsize*2;

for iiEdge = 1:length(edgeLengths)
    if edgeLengths(iiEdge)>lengthCutoff
       rejectBeads(edgesList(iiEdge,1)) = 1; 
       rejectBeads(edgesList(iiEdge,2)) = 1;
    end
end

%%
for iiBead = 1:size(beadStruct.psfExample,2)
%     if not(rejectBeads(iiBead))
        tempPSF = beadStruct.psfExample(iiBead).psfCutOut;
    surfacePSFcnts = surfacePSFcnts+not(isnan(tempPSF));
    tempPSF(isnan(tempPSF)) = 0;
    surfacePSF = surfacePSF+tempPSF;
%     nBeads = nBeads +1;
%     end
end

end

%%
surfacePSFnorm = surfacePSF./surfacePSFcnts;
surfacePSFnorm = surfacePSFnorm./max(surfacePSFnorm(:));
SliceBrowser(surfacePSFnorm);

% surfPsf2 = (surfacePSFnorm+flipdim(surfacePSFnorm,1)+...
%     flipdim(surfacePSFnorm,2)+flipdim(flipdim(surfacePSFnorm,1),2))/4;
% % SliceBrowser(surfacePSFnorm-surfPsf2);
% surfPsf3 = bpass3(surfacePSFnorm,0.5,[],-1);
% surfPsf3 = surfPsf3./max(surfPsf3(:));
% % weight the mixing of the two psfs by this
% % surfPsf4 = exp(-abs(surfPsf2)./0.03).^3;
% % surfPsf5 = surfPsf4.*surfPsf3+(1-surfPsf4).*surfPsf2;
% % 
% SliceBrowser(surfPsf5-surfPsf3);
%%
x = linspace(0,1,100);
y = exp(-x/.03);
plot(x,y);
%% (2b) group PSFs into sections, get set up
zPSFsectionCenter = -5:0.3:1;
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
    if size(psfExample,2)<20
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

