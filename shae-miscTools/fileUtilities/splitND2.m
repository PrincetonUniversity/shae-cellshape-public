function splitND2();
%%
[fileName,pathName] = uigetfile('*.nd2');
%%
fileData = bfopen(fullfile(pathName,fileName));
%% parse first label to determine number of images
nSeries = size(fileData,1);

for jSeries = 1:nSeries
   % each series may be stored different
   tempString = fileData{jSeries}{1,2};
    clear concatImage colorLabels;
    [iSeries,nSeries,iPlane,nPlanes,iColor,nColors,iTime,nTimes,stackString] =...
        parseND2label(tempString);
    for jPlane = 1:nPlanes
        tempString = fileData{jSeries}{jPlane,2};
        [iSeries,nSeries,iPlane,nPlanes,iColor,nColors,iTime,nTimes,stackString] =...
            parseND2label(tempString);
        concatImage(:,:,iTime,iColor) = fileData{jSeries}{jPlane,1};
        if iTime ==1
            colorLabels{iColor}=[stackString,'.tif'];
        end
    end
    
    for iColor = 1:length(colorLabels)
       tiffwrite(fullfile(pathName,colorLabels{iColor}), concatImage(:,:,:,iColor));
    end
end