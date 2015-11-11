function info = parseTags(tifFilePath)
% This function attempts to open a text file with the same file stem as a
% .tif file and parse the information within it. This function needs to be
% updated whenever the labview code generating the text file changes.
% 
% The text file should have the information structured as
% fieldname:fieldvalue NEWLINE
%
% Some fields are converted from their string value into the proper
% datatypes. These fields are
% timestamp -> matlab datenum
% exposureTime -> double
% posX, posY, posZ -> double
%
% 20111206 BPB first version

%% open .txt file with metaData and parse it into tags
FIDa = fopen(strrep(tifFilePath,'.tif','.txt'));

info = struct;
if FIDa<3
    return
end

while not(feof(FIDa))
    infoTemp = fgetl(FIDa);
    fieldnameBreak = strfind(infoTemp,':');
    if not(isempty(fieldnameBreak))
        fieldnameBreak = fieldnameBreak(1);
        fieldnameString = infoTemp(1:fieldnameBreak-1);
        fieldnameString = textscan(fieldnameString,'%s','Delimiter',' \b\t=()');
        fieldnameString2 = [];
        for iWord=1:size(fieldnameString{1},1);
            fieldnameString2 = [fieldnameString2,fieldnameString{1}{iWord}];
        end
        
        info = setfield(info,fieldnameString2,infoTemp(fieldnameBreak+1:end));
    end
           
        
end
fclose(FIDa);

%% clean up fields that are not actually strings
% date and time
if isfield(info,'timestamp')
    info.timestamp = datenum(info.timestamp,'yyyy-mm-dd_HH:MM:SS.FFF');
end

% positions
if isfield(info,'posX')
    info.posX = str2double(info.posX);
end
if isfield(info,'posY')
    info.posY = str2double(info.posY);
end
if isfield(info,'posZ')
    info.posZ = str2double(info.posZ);
end

% exposure time
if isfield(info,'exposureTime')
    info.exposureTime = str2double(info.exposureTime);
end

% cycle time
if isfield(info,'KineticCycleTime')
    info.KineticCycleTime = str2double(info.KineticCycleTime);
end

% FRAP/FLIP fields
if isfield(info,'pulseDelay');
    info.pulseDelay = str2double(info.pulseDelay);
end
if isfield(info,'pulseDuration');
    info.pulseDuration = str2double(info.pulseDuration);
end