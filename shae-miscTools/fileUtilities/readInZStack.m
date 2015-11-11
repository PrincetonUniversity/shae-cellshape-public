function [zStackOut,varargout] = readInZStack(folderName,param)
% original version BPB 2013
if isappdata(0,'readInZStackDefaultFolder')
    defaultFolder = getappdata(0,'readInZStackDefaultFolder');
else
    defaultFolder = pwd;
end

% ask for folder to read in
if nargin<1
    [folderName] = uigetdir(defaultFolder,'Please select the folder to read in');
end

if nargin<2
    param = struct();
end

if not(isfield(param,'isDisplayStatusBar'))
    param.isDisplayStatusBar = true;
end

if not(isfield(param,'complicatedStatusBar'))
   param.complicatedStatusBar = 'none'; 
end

% store parent directory for fewer clicks
defaultFolder = fileparts(folderName); % selected folder
defaultFolder = fileparts(folderName); % parent directory
setappdata(0,'readInZStackDefaultFolder',defaultFolder);
%% prepare to read in stack
fileList = dir(fullfile(folderName,'*.tif*'));

% read in each file
zStackOut.imageInfo = {[]};
imageSize = size(imread(fullfile(folderName,fileList(1).name)));
zStackOut.image = nan([imageSize,length(fileList)]);
zStackOut.zPos = [];
%% read in stack
if param.isDisplayStatusBar
    progressbar(0);
end

if strcmp(param.complicatedStatusBar,'none')
else
    % the complicated bar is passed in as pairs of numbers, each pair is
    % the appropriate 0 to 1 scaling for the subsegment below
    
    % for example
%     complicatedStatusBar = [[0,0.3];[0,1]]; % two level bar, as the second bar
%     completes, the first bar goes from 0 to 0.3
%
%     complicatedStatusBar = [[0,0.3];[0,0.25];[0,1]] % three level bar
%      third bar goes 0 to 1
%      second bar goes 0 to 0.25
%      first bar goes 0 to 0.25*0.3=0.075
%
%     complicatedStatusBar = [[0,0.3];[0.5,0.6];[0.75,1.0];[0,1]] % four level bar
%      fourth bar goes 0 to 1
%      third bar goes 0.75 to 1
%      second bar goes 0.75*(0.6-0.5)+0.5=0.575 to 0.6
%      first bar goes for 0.3*(0.575) to 0.3*(0.6)
%     

    currentPosition = param.complicatedStatusBar(:,1);
    totalStep = diff(param.complicatedStatusBar,1,2);
    nDimBar = size(currentPosition,1);
    % figure out how much it needs to increased by
    progBarStepSize(nDimBar) = (1/length(fileList));
    for iiDimBar = nDimBar-1:-1:1
        % cludgy notation, if the bar below this one doesn't change, ignore it and
        % check the next
        if totalStep(iiDimBar+1) == 0;
            progBarStepSize(iiDimBar) = progBarStepSize(iiDimBar+2).*totalStep(iiDimBar);           
            
        else
            progBarStepSize(iiDimBar) = progBarStepSize(iiDimBar+1).*totalStep(iiDimBar);
        end
    end

    
    buildProgressBarString(currentPosition)

end


for iFile = 1:length(fileList);
    % note: this construction can return the original array without any
    % filtering
    zStackOut.image(:,:,iFile) = bpass(double(imread(...
        fullfile(folderName,fileList(iFile).name))),...
        0,0);
    
    zStackOut.imageInfo{iFile} = parseTags(  ...
        fullfile(folderName,fileList(iFile).name));
    if not(isempty(fieldnames(zStackOut.imageInfo{iFile})))
        zStackOut.zPos(iFile) = zStackOut.imageInfo{iFile}.posZ;
    end
    
    if param.isDisplayStatusBar
        progressbar(iFile/length(fileList));
    end
    
    if strcmp(param.complicatedStatusBar,'none')
    else
            currentPosition = currentPosition+progBarStepSize';
            buildProgressBarString(currentPosition);
    end
end

if isempty(zStackOut.zPos)
    sliceSize = inputdlg('Spacing between planes in um?');
    zStackOut.zPos = (1:size(zStackOut.image,3)).*str2double(sliceSize);
end


