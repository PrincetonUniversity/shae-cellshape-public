function readInCsv_signalSingleCells();
% for datasets with large numbers of cells, this will read back in the
% output csv files for generating plots with different datasets
csvFileList = uipickfiles();
%%
summaryCellArray = cell(1,4);
for iiCsvFile = 1:length(csvFileList);
    T = readtable(csvFileList{iiCsvFile},'Delimiter',',');
    %%
    % final summary table is a cell array with nFolders x 4 columns
    
    % column 1 is the full path to the source folder
    summaryCellArray{iiCsvFile,1} = 'reread from csv';
    
    
    %%
    try
        % column 2 is the channel that was brighter
        summaryCellArray{iiCsvFile,2} = T.maxChannelIDX;
        
        % column 3 is the list of measured intensities (variable number of
        % channels)
        channelIntensities = table2array(T(:,3:end));
    catch ME1
            
        warning('signalSingleCells:readInCsv:legacyFormat',...
            'Legacy csv detected. Assuming first column with filename is missing.');
        % column 2 is the channel that was brighter
        summaryCellArray{iiCsvFile,2} = T.Var1;
        
        % column 3 is the list of measured intensities (variable number of
        % channels)
        channelIntensities = table2array(T(:,2:end));
        
    end
    summaryCellArray{iiCsvFile,3} = channelIntensities;
    
    % column 4 is the name of the source file (or the output file)
    [saveFolderName,saveFileName,saveFileExt] = fileparts(csvFileList{iiCsvFile});
    summaryCellArray{iiCsvFile,4} = [saveFileName,saveFileExt];
    
    % remove rows that have a negative value
    removeRow = any(channelIntensities<0,2);
    if any(removeRow)
       warning('signalSingleCells:readInCsv:negativeValues',['Negative intensity values. Removing ',...
          num2str(nnz(removeRow)), ' cell(s) from analysis']);
      channelIntensities(removeRow,:) = [];
      summaryCellArray{iiCsvFile,3} = channelIntensities;
      tempMaxChannelIDX = summaryCellArray{iiCsvFile,2};
      tempMaxChannelIDX(removeRow) = [];
       summaryCellArray{iiCsvFile,2} =      tempMaxChannelIDX;
    end
    
    %%
    if iiCsvFile == 1;
        nSplits = size(channelIntensities,2);
        maxSignalTotal = max(channelIntensities);
    else
        nChannelsII = size(channelIntensities,2);
        if nChannelsII~=nSplits
            error('signalSingleCells:readInCsv:mismatchSize',...
                'Mismatch in number of channels. Please select csv files with same number of channels or contact the developer.');
        else
            maxIntensityTemp = max(channelIntensities);
            maxSignalTotal = max(maxSignalTotal,maxIntensityTemp);
        end
    end
    
end
%%
displayPlots_signalSingleCells(summaryCellArray,nSplits,maxSignalTotal)


