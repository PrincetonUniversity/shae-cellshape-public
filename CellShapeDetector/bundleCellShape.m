function bundleCellShape(folderPath);
%%
if nargin<1
    folderPath = pwd;
end
%disp(folderPath);
%% grab everyone in the main folder
cellShapePath = which('Cell_shape_detector3dConvFinalFolder.m');
[cellShapePath] = fileparts(cellShapePath);
cellShapeList = [];

% find the directory just above CellShapeDetector
[csParent,cellShapeName] = fileparts(cellShapePath);
[csGParent,csParent] = fileparts(csParent);

% folderPath

for ext = [{'*.m'},{'*.mat'},{'*.fig'},{'*.p'},{'*.mex'}];
    cellShapeListTemp = dir([cellShapePath,filesep,ext{1}]);
    for jj = 1:numel(cellShapeListTemp)
        if strcmp(cellShapeListTemp(jj).name(1:2),'._');
        else
            tempName = {[cellShapePath,filesep,cellShapeListTemp(jj).name]};
            cellShapeList = cat(1,cellShapeList,tempName);
        end
    end
end
% grab everyone thought to be a dependent
cellShapeList2 = (mydepfun('Cell_shape_detector3dConvFinalFolder',1));


% smush these two lists together
cellShapeList = unique(cat(1,cellShapeList,cellShapeList2));

% see if any are not in the correct place
for ll = numel(cellShapeList):-1:1
    tempName = cellShapeList{ll};
    [tempParent,tempName,tempExt] = fileparts(tempName);
    [tempGParent,tempParent] = fileparts(tempParent);
    [tempGGParent,tempGParent] = fileparts(tempGParent);
    if not(strcmp(tempGParent,csParent))
        warning('bundleCellShape:extraFiles',['Extra file found: ',...
            cellShapeList{ll},'. Not included in file listing.']);
        cellShapeList{ll} = [];
    else
        cellShapeList{ll} = [tempParent,filesep,tempName,tempExt];
    end
end
%
%% write out list of files
fidA  = fopen(fullfile(folderPath,'cellShape.txt'),'w');
for ii = 1:numel(cellShapeList)-1;
   fwrite(fidA,[cellShapeList{ii}]);
   fprintf(fidA,'\r\n');
end
% don't add a new line to the last one
   fwrite(fidA,[cellShapeList{ii+1}]);
% close the file
fclose(fidA);
%%
% zip up files
% zip(fullfile(folderPath,'cellShape.zip'),cellShapeList);
% fprintf(1,'Created %s with %d entries\n',fullfile(folderPath,'cellShape.zip'),numel(cellShapeList));

% exportToZip('Cell_shape_detector3dConvFinalFolder.m',fullfile(folderPath,'cellShape.zip'));