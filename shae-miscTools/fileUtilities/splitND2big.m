function splitND2big(filepath)
%%
if nargin<1
    [fileName,pathName] = uigetfile('*.nd2');
else
    % parse this back into the appropriate fields
    [pathName,fileName,fileExt] = fileparts(filepath);
    pathName = [pathName,fileExt];
end
%%
% this function is similar to splitND2 except it assumes that the image
% will not be able to be held in memory

mPlanesPerImage = inputdlg('How many planes to save in an image? -1 for all');

mPlanesPerImage = str2double(mPlanesPerImage);
if isnan(mPlanesPerImage)
    mPlanesPerImage = 0;
end

if mPlanesPerImage>0
    nOverlap =  inputdlg('How many planes to overlap between image files?');
    nOverlap = str2double(nOverlap);
    if isnan(mPlanesPerImage)
        nOverlap = 0;
    end
end
%%
reader = bfGetReader(fullfile(pathName,fileName));
numSeries = reader.getSeriesCount;
id = fullfile(pathName,fileName);

% fileData = ;
%% parse first label to determine number of images

for jSeries = 1:numSeries
    % Set the current series index
    reader.setSeries(jSeries-1);
    % Find out how many images this series has
    numImages = reader.getImageCount;
    metadataList = reader.getMetadata();
    metadataStore = reader.getMetadataStore();
        
    % from bfopen
    % build an informative title for our figure
    imageIndexI = 1;
    label = id;
    if numSeries > 1
        qs = int2str(jSeries);
        label = [label, '; series ', qs, '/', int2str(numSeries)];
    end
    
    qi = int2str(imageIndexI);
    label = [label, '; plane ', qi, '/', int2str(numImages)];
    if reader.isOrderCertain()
        lz = 'Z';
        lc = 'C';
        lt = 'T';
    else
        lz = 'Z?';
        lc = 'C?';
        lt = 'T?';
    end
    
    zct = reader.getZCTCoords(imageIndexI - 1);
    sizeZ = reader.getSizeZ();
    if sizeZ > 1
        warning('splitND2big:zStacks','This function does not yet support multiple zplanes');
        qz = int2str(zct(1) + 1);
        label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
    end
    sizeC = reader.getSizeC();
    if sizeC > 1
        qc = int2str(zct(2) + 1);
        label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
    else
        % sometimes there won't be a color channel
        qc = int2str(0);
    end
    sizeT = reader.getSizeT();
    if sizeT > 1
        qt = int2str(zct(3) + 1);
        label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
    end
    
    
    % in lieu of parseND2label
    iColor = str2double(qc);
    stackString = [fileName,'_series',num2str(jSeries),'_C=',num2str(iColor)];
    stackString = strrep(stackString,'.nd2','');
    nColors = sizeC;
    nTimes = sizeT;
    
%     % each series may be stored different
%     clear concatImage colorLabels;
%     [iSeries,numSeries,iPlane,nPlanes,iColor,nColors,iTime,nTimes,stackString] =...
%         parseND2label(label);
%     % BPB finds this humorous, the parseND2label is merely a way of undoing
%     % the bfopen generation of a label for the image
%     
    
    % divide into mPlanesPerImage stacks
    kTime = 1;
    % start with single elements, add on as the image is bigger
    timeStart = 1;
    timeEnd = 1;
    kImage = 1;
    % list of planes in each image
    while kTime<=nTimes
        kTime = kTime+1;
        timeEnd(kImage) = kTime-1;
        % then the overlap case
        if (timeEnd(kImage)-timeStart(kImage))== mPlanesPerImage;
           kImage = kImage+1;
           kTime = kTime+1-nOverlap;
           timeStart(kImage) = kTime;
        end
    end
    
    stackString(strfind(stackString,'C='):end) = [];
    for iTime = 1:length(timeStart)
        % first grab the names of the color channnels
        for iColor = 1:nColors
            tiffFileName{iColor,iTime} = fullfile(pathName,...
                [stackString,num2str(timeStart(iTime)),'-',num2str(timeEnd(iTime)),...
                '_','C=',num2str(iColor-1),'.tif']);
            if exist(tiffFileName{iColor,iTime})
               warning('splitND2big:deletingFile',...
                   ['Deleting file: ', tiffFileName{iColor,iTime}]);
               delete(tiffFileName{iColor,iTime});
            end
        end
    end
    
    % indices for time are 0 based?
    timeStart = timeStart - 1;
    timeEnd = timeEnd - 1;
    
    % loop over each subimage
    for kImage = 1:length(timeEnd)
        % loop over each image
        for imageIndexI = 1:(numImages);
            % build a label
            label = id;
            if numSeries > 1
                qs = int2str(jSeries);
                label = [label, '; series ', qs, '/', int2str(numSeries)];
            end
            
            qi = int2str(imageIndexI);
            label = [label, '; plane ', qi, '/', int2str(numImages)];
            if reader.isOrderCertain()
                lz = 'Z';
                lc = 'C';
                lt = 'T';
            else
                lz = 'Z?';
                lc = 'C?';
                lt = 'T?';
            end
            zct = reader.getZCTCoords(imageIndexI - 1);
            sizeZ = reader.getSizeZ();
            if sizeZ > 1
                warning('splitND2big:zStacks','This function does not yet support multiple zplanes');
                qz = int2str(zct(1) + 1);
                label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
            end
            sizeC = reader.getSizeC();
            if sizeC > 1
                qc = int2str(zct(2) + 1);
                label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
            end
            sizeT = reader.getSizeT();
            if sizeT > 1
                qt = int2str(zct(3) + 1);
                label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
            end
            
            % if this is an image in the proper stack, drop it into place
            isInSubimage =  and(zct(3)>=timeStart(kImage),zct(3)<=timeEnd(kImage));
            
            if isInSubimage
                tempImage = bfGetPlane(reader,imageIndexI);
                imwrite(tempImage,tiffFileName{zct(2)+1,kImage},'tif','WriteMode','append')
                if mod(imageIndexI,10)==0
                    fprintf(1,'.');
                end
            end
            
            
        end
    end
    
end
    
disp('done');

    %%%
    %%
    %%%
% % 
% %     counter = 0;
% %     kImage = 1;
% %     jPlane = 1;
% %     
% %     
% %     
% %     while jPlane<nPlanes
% %         
% %         % increment counters
% %         counter = counter+1;
% %         %
% %         if counter == mPlanesPerImage;
% %             %
% %             kImage = kImage + 1;
% %             counter = 1;
% %             
% %             for kColor = 1:nColors
% %                 % new file name
% %                 tiffFileNameNew{kColor} = strrep(tiffFileName{kColor},'C=',[num2str(1+((jPlane-1)/nColors)),'C=']);
% %                 
% % %                 % check if this new file already exists
% % %                 if exist(tiffFileNameNew{kColor},'file');
% % %                     warning('splitND2Big:removePreviousFile',...
% % %                         ['Overwriting previous file', tiffFileNameNew{kColor}]);
% % %                     delete(tiffFileNameNew{kColor});
% % %                 end
% %                 movefile(tiffFileName,tiffFileNameNew);
% %             end
% %             
% %             % move back a few frames if overlap is desired
% %             jPlane = jPlane - nColors*nOverlap;
% %         end
% %         
% %         tempString = fileData{jSeries}{jPlane,2};
% %         [iSeries,nSeries,iPlane,nPlanes,iColor,nColors,iTime,nTimes,stackString] =...
% %             parseND2label(tempString);
% %         %         concatImage(:,:,iTime,iColor) = fileData{jSeries}{jPlane,1};
% %         tempImage = fileData{jSeries}{jPlane,1};
% %         % need to write out as a file, with append mode
% %         
% %         if subImageIndex == 1;
% %             tiffFileName{iColor} = fullfile(pathName,...
% %                 [strrep(stackString,'C=',[num2str(planeIndex),'C=']),'.tif']);
% %         end
% %                 
% % %         if counter == 1;
% % 
% % %             %
% % %             for kColor = 1:nColors
% % %                 % new file name
% % %                 tiffFileNameNew{kColor} = strrep(tiffFileName{kColor},'.tif',[num2str(jPlane/nColors)),'.tif']);
% % %                 
% % %                 % check if this new file already exists
% % %                 if exist(tiffFileNameNew{kColor},'file');
% % %                     warning('splitND2Big:removePreviousFile',...
% % %                         ['Overwriting previous file', tiffFileNameNew{kColor}]);
% % %                     delete(tiffFileNameNew{kColor});
% % %                 end
% % %                 movefile(tiffFileName,tiffFileNameNew);
% % %             end
% % %             
% % %             % move back a few frames if overlap is desired
% % %             jPlane = jPlane - nColors*nOverlap;
% % %         end
% %         
% %         imwrite(tempImage,tiffFileName{iColor},'tif','WriteMode','append');
% %         %
% %         jPlane = jPlane+1;
% %         
% %     end
% %     
% %     
% % %     if counter == 1;
% % %         tiffFileNameNew = strrep(tiffFileName,'.tif',[num2str(1+((jPlane-1)/nColors)),'.tif']);
% % %         if exist(tiffFileNameNew,'file')
% % %             warning('splitND2Big:removePreviousFile',...
% % %                 ['Overwriting previous file', tiffFileNameNew]);
% % %             pause();
% % %             delete(tiffFileNameNew);
% % %         end
% % %         movefile(tiffFileName,tiffFileNameNew);
% % %     end
% %     
% %     %     for iColor = 1:length(colorLabels)
% %     %        tiffwrite(fullfile(pathName,colorLabels{iColor}), concatImage(:,:,:,iColor));
% %     %     end
% % end
