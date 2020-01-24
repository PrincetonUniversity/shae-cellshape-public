function measureSignalSingleCroppedCells()

%%
folderList = uipickfiles();

saveFolder = uigetdir(folderList{1},'Please select a folder to store results in.');
%%

% initialization
nSplits = 2; % how many total channels are there in each image?

maxSignalTotal = [0,0]; % continuously increase the maxSignal as new images are read


% loop over all folders
progressbar('Folders','Files');
nFolders = length(folderList);
for iiFolder = 1:nFolders;
    % fraction of the folders finished
    folderFraction = (iiFolder-1)./nFolders+100*eps;
    progressbar(folderFraction,0);
    
    % parse folder names as descriptions
    [parent,folder] = fileparts(folderList{iiFolder});
    [gparent,parent] = fileparts(parent);
    [ggparent,gparent] = fileparts(gparent);
    
    descriptiveString = [gparent,'_',parent,'-',folder,'_',transcodeDate(now),'.csv'];
    
    % open a csv file
    fidA = fopen(fullfile(saveFolder,descriptiveString),'w+');
    fprintf(fidA,['Filename,','maxChannelIDX,']);
    for iiChannel = 1:(nSplits-1)
        fprintf(fidA,['Signal_Channel_',num2str(iiChannel),',']);
    end
    fprintf(fidA,['Signal_Channel_',num2str(nSplits),'\r']);
    
    % generate a list of files to process
    fileList = dir([folderList{iiFolder},filesep,'*.tif']);
    
    for iiFile = 1:length(fileList)
        % update progressbar
        fileFraction = (iiFile-1)/length(fileList)+100*eps;
        folderFraction2 = folderFraction + fileFraction.*(1/nFolders);
        progressbar(folderFraction2,fileFraction);
        fullFilePath = fullfile(folderList{iiFolder},fileList(iiFile).name);
        fprintf(fidA,[fileList(iiFile).name,',']);

        % read in and process image
        tiffArray = readImage(fullFilePath);
        [tempMaxIDX,tempSignal] = processTiffArray(tiffArray,nSplits);
        if any(tempSignal<0)
           warning('signalSingleCell:measureSignalSingleCropped:negativeIntensity',...
               ['Negative intensity detected in ', fullFilePath ,'.']);
        end
        signal(iiFile,:) = tempSignal;
        maxSigChannelIDX(iiFile) = tempMaxIDX;

        
        fprintf(fidA,[num2str(maxSigChannelIDX(iiChannel)),',']);
        
        % write out to csv
        for iiChannel = 1:(nSplits-1)
            fprintf(fidA,[num2str(signal(iiFile,iiChannel)),',']);
        end
        
        fprintf(fidA,[num2str(signal(iiFile,nSplits)),'\r']);
        
    end
    
    
    
    % close file for reading later
    fclose(fidA);
    
    % summarize into a cell array
    summaryCellArray{iiFolder,1} = folderList{iiFolder};
    summaryCellArray{iiFolder,2} = maxSigChannelIDX;
    summaryCellArray{iiFolder,3} = signal;
    summaryCellArray{iiFolder,4} = descriptiveString;
    maxSignalTotal = max(maxSignalTotal,max(signal));
end
progressbar(1);
displayPlots_signalSingleCells(summaryCellArray,nSplits,maxSignalTotal);


function tiffArray = readImage(fullFilePath);
% testing on 20190614 R2018b,this style of reading in is about 35% faster
% than tiffread for the types of images used (~100x100x50 voxels)

% file storage information
fileInfo = imfinfo(fullFilePath);
imageWidth = fileInfo(1).Width;
imageHeight = fileInfo(1).Height;
nPlanes = numel(fileInfo);

if fileInfo(1).BitDepth == 32
    % Single precision floating point tiff, preallocate, otherwise, don't
    % preallocate
    
    % testing on 20190614 R2018b, preallocated with nans is ~3% faster than letting
    % matlab do it on the fly
    tiffArray = nan(imageHeight,imageWidth,nPlanes,'single');
else
    
end

for iiPlane = 1:size(fileInfo);
    tiffArray(:,:,iiPlane) = imread(fullFilePath,'index',iiPlane,'info',fileInfo);
end

function [maxSigChannelIDX,signal] = processTiffArray(tiffArray,nSplits);
reshapeImage = reshape(tiffArray,size(tiffArray,1),size(tiffArray,2),[],nSplits);
sumProj = (sum(reshapeImage,3));

isDisplayFigures = false;
if isDisplayFigures
    figure(4);
    clf;
    for iiChannel = 1:nSplits
        subplot(nSplits,1,iiChannel)
        imshow(sumProj(:,:,1,iiChannel),[]);
        title(['channel ',num2str(iiChannel),': ',num2str(sum(sum(sumProj(:,:,1,iiChannel))))]);
    end
end

% segment based on whichever channel has more total signal
totalSignal = squeeze(sum(sum(sumProj,1),2));
[maxSig,maxSigChannelIDX] = max(totalSignal);

% segment image
templateImage = sumProj(:,:,1,maxSigChannelIDX);
templateImage = templateImage-min(templateImage(:));
templateImage = templateImage./max(templateImage(:));
%    [counts,x] = imhist(templateImage,512);
%    threshLevel = otsuthresh(counts);
threshImage = imbinarize(templateImage);
threshImage = bwmorph(threshImage,'clean');
threshImage = bwmorph(threshImage,'dilate');

% calculate the signal in each channel in these regions
for iiChannel = 1:nSplits
    signal(iiChannel) = sum(sum(threshImage.*sumProj(:,:,1,iiChannel),1),2)./...
   nnz(sumProj(:,:,1,iiChannel));
end
