%%This is the main operating code for cell shape detector. It fits the
%3d stack of a cell with the contour that best describes the cell surface
%after convolution with PSF.
function Cell_shape_detector3dConvTriFolder(folderName,folderIdx,fileGroupSize)

%% to turn off specific warnings
% [warnmsg,warnid] = lastwarn()
% warning('off',warnid);


%% exit codes
% 0 -
% 1 -
% 2 - no settings files found
% 3 - folderName is not found
tic
%% process the input arguments to determine the mode of running
% file group size is how many files to group together in one run of the
% main code
if nargin<3
    fileGroupSize = 1;
end

% sometimes these variables are passed in as strings, sometimes as numbers
if ischar(fileGroupSize)
    fileGroupSize = str2double(fileGroupSize);
end

if ischar(folderIdx)
    folderIdx = str2double(folderIdx);
end

%% usage type 2, folderIdx is a nan and the input folder name is actually a file
disp(folderName);
disp(folderIdx);
% check to see if this is a file or a directory
if exist(folderName,'dir') 
    % if it is a directory, probably usage type 1
    usageType = 1;
elseif exist(folderName,'file');
    % this a file, usage type two. 
    usageType = 2;
    [folderName,fileName,ext] = fileparts(folderName);
    fileName = [fileName,ext];
    
    % at a later point, one will need to find which fileIDX this fileName
    % corresponds with
else
    % this is not a folder nor a file
    disp('Exiting code 3: neither a folder nor file.')
    disp(['folderName: ', folderName]);
    disp(['pwd: ',pwd]);
%     exit(3);
end


%
if fileGroupSize>1
    if numel(folderIdx)>1
        warning('cellShape:invalidParams','Suppressing extra folderIdx''s and using fileGroupSize');
    end
    folderIdx = folderIdx(1):folderIdx(1)+fileGroupSize-1;
    
end

disp(['Processing : ',folderName,'. Files : ',num2str(folderIdx),'.']);


try
    %if tri type settings txt exist, use it, otherwise use whatever is there
    if exist([folderName filesep 'cell_shape_settings_tri.txt'],'file')
        [image_para,flag,cline_para] = paramsText(fullfile(folderName,'cell_shape_settings_tri.txt'));
        
    else
        [image_para,flag,cline_para] = paramsText(fullfile(folderName,'cell_shape_settings.txt'));
    end
catch ME
    save('errorCaught.mat');
%     disp(fullfile(folderName,'cell_shape_settings_tri.txt'));
%     disp(pwd)
%     keyboard
%    exit(2); 
end

disp(['Processing params: ',num2str(toc),' seconds.']);
%%


% string to append on end of output
ptail='TRI.mat';

% get a list of files in this directory
fileList = dir([folderName,filesep,'*.tif']);

% remove annoying ._
fileList=fileList(cellfun(@(x) isempty(strfind(x,'._')),{fileList.name}));

switch usageType
    case 1 
        % don't do anything, folderIdx is already a number and the 
        % identity in the list
    case 2
       folderIdx = find(strcmp(fileList,fileName));
end

%% loop over files ------

for kkFile = 1:numel(folderIdx);
    ticA = tic;
    p1 = folderName;
    
    clearvars -except image_para flag cline_para ptail folderName folderIdx kkFile fileList p1 ticA
    stack_name = fileList(folderIdx(kkFile)).name;
    [~,stack_name,~] = fileparts(stack_name);
    disp(['Starting on ',stack_name,'.']);
    try
        
        % full file path to name of image
        %p=strcat(p1,filesep,strrep(stack_name,'.tif',''));
        p=strcat(p1,filesep,(stack_name)); %removing extra copies of '.tif' can cause issues
        
        % ---- Parameters that are not used within structures
        pixel_interp = flag.pixel_interp;
        intensity_correction = flag.intensity_correction;
        window = flag.window;
        end_angles = flag.end_angles;
        pixel_interp = flag.pixel_interp;
        Fflag2 = flag.Fflag2;
        % ----
        
        %% Obtaining size of im
        [row col] = size(imread([p,'.tif'], 'tif', 1));
        
        % If there are multiple time points and there is a fluorescent channel to
        % be mapped, get a folder ready for them
        %if image_para.Fstack_t_size>1 && flag.F1==1
        %    mkdir(p)
        %end
        
        %% load raw hyper stack
        
        %hyper_stackx = zeros(row, col);
        
        if ~flag.F1
            image_para.stack_z_size=min(length(imfinfo([p,'.tif'])),...
                image_para.stack_z_size);
        end
        
        
        %% find and load metadata
        %%% ----- filter shift patch
        
        % find the last dash in the filename which separates images into individual cells
        % this dash has ascii code 45
        try
            dashPos = find(double(stack_name)==45,1,'last');
            metaFileName = stack_name;
            metaFileName(dashPos:end) = [];
            metaFileName = [metaFileName,'-Metadata.mat'];
            load(metaFileName);
            isMetaData = 1;
        catch ME
            isMetaData = 0;
        end
        
        %%% ----- end filter shift patch
        %% load raw hyper stack
        
        hyper_stack=zeros(row,col,image_para.stack_z_size);
        
        for k = 1:image_para.stack_z_size
            % determine which plane should be read in
            if flag.reverse==1 && flag.F1==1
                planeIndex = image_para.Fstack_t_size*image_para.Fstack_z_size+k;
            else
                planeIndex = k;
            end
            % read in the plane and store it
            hyper_stack(:,:,k) = double(imread([p,'.tif'], 'tif',...
                planeIndex)...
                +intensity_correction);
            %  hyper_stack(:,:,k)=wiener2(hyper_stack(:,:,k));
            %             if ~isempty(find(hyper_stack(:,:,k)<0, 1)) && k==1
            %                 H1=hyper_stack(:,:,k);
            %                 Hmask=(H1==0);
            %                 Hmask(smooth2a(Hmask,3,3)>0)=1;
            %
            %             elseif ~isempty(find(hyper_stack(:,:,k)<0, 1))
            %                 H1=hyper_stack(:,:,k);
            %                 Hmin=min(H1(~Hmask));
            %                 H1(Hmask)=Hmin;
            %                 hyper_stack(:,:,k)=H1;
            %             end
            
            
        end
        
        %%% ----- filter shift patch
        if isMetaData
            % background subtract, shiftImgFilter, background add
            backgndTmp = quantile(hyper_stack(:),0.005);
            hyper_stack = hyper_stack-backgndTmp;
            for k2 = 1:image_para.stack_z_size
                if flag.reverse==1 && flag.F1==1
                    planeIndex = image_para.Fstack_t_size*image_para.Fstack_z_size+k2;
                else
                    planeIndex = k2;
                end
                
                try
                filterString = totalImageInfo{planeIndex}.filterID;
                catch ME % sometimes there is partial metadata but not all of it
                    filterString = [nan,nan];
                end
                filterIdx = str2double(filterString(end-1));
                if isnan(filterIdx)
                    filterIdx = 3; % no shifting of image
                end
                
                hyper_stack(:,:,k2) = shiftImgFilter(hyper_stack(:,:,k2),filterIdx);
                
            end
            hyper_stack = hyper_stack+backgndTmp;
        end
        %%% ----- end filter shift patch
        
        hyper_stack=image_resize(hyper_stack,size(hyper_stack,1),...
            size(hyper_stack,2),round(image_para.ZX_factor*size(hyper_stack,3)));
        image_para.imsize=size(hyper_stack);
        pedMask=false(image_para.imsize);
        pedMask(1:3,:,:)=true;
        pedMask(:,1:3,:)=true;
        pedMask(end-2:end,:,:)=true;
        pedMask(:,end-2:end,:)=true;
        pedMask(hyper_stack==0)=false;
        
        pedistal=median(hyper_stack(pedMask));
        hyper_stack(hyper_stack<pedistal)=pedistal;
        hyper_stack=normalizeRange(hyper_stack);
        
        
        
        
        %% Paramaterize centerline by distance instead of coordinate
%% Paramaterize centerline by distance instead of coordinate
%  z_level=.3;
 
% [xyzs,zbw,zsq]=cellCenterline3d(hyper_stack,z_level,window,...
%     image_para.imsize(3),flag);
% 
% %% Calculate tangent vectors of points and create new basis for detecting membrane
% [x,y,z] = ndgrid(-13:13);
% CLim=coord2image3d(xyzs(:,1),xyzs(:,2),xyzs(:,3),image_para.imsize,...
%     1,1);
% CLim=imdilate(CLim,strel(sqrt(x.^2 + y.^2 + z.^2) <=cline_para.radius));
if ~isfield(cline_para,'upperLimit')
    cline_para.upperLimit=2;
    cline_para.lowerLimit=.5;
end


% bandpass in 3D
filt3DStack = bpass3_jn(hyper_stack,.25,[7,7,30]);
filt3DStack=normalizeRange(filt3DStack);
for ii=1:size(filt3DStack,3)
    temp=bwmorph(...
    normalizeRange(filt3DStack(:,:,ii))>flag.z_level,'clean');
temp=temp.*(filt3DStack(:,:,ii)>flag.z_level/2);

filt3DMask(:,:,ii)=temp;
end
%%
% threshold
%filt3DMask = (filt3DStack./max(filt3DStack(:)))>z_level;
filt3DMask= imclearborder(filt3DMask,6);

filt3DMask(1:2,:,:)=0;
filt3DMask(end-1:end,:,:)=0;
filt3DMask(:,1:2,:)=0;
filt3DMask(:,end-1:end,:)=0;
filt3DMask(:,:,1:2)=0;
filt3DMask(:,:,end-1:end)=0;

CC=bwconncomp(filt3DMask,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = numPixels==max(numPixels);
for Iidx=1:length(idx)
    if ~idx(Iidx)
filt3DMask(CC.PixelIdxList{Iidx}) = 0;
    end
end

[x,y,z] = ndgrid(-9:9);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=9);
se2 = strel(sqrt(x.^2 + y.^2 + z.^2) <=15);


filt3DMask= imclearborder(filt3DMask);
filt3DMask= imdilate(filt3DMask,se);
for runs=1:2
for imSlice=1:size(filt3DMask,1);
    filt3DMask(imSlice,:,:)=imfill(squeeze(filt3DMask(imSlice,:,:)),'holes');
end    
for imSlice=1:size(filt3DMask,3);
    filt3DMask(:,:,imSlice)=imfill(filt3DMask(:,:,imSlice),'holes');
end
for imSlice=1:size(filt3DMask,2);
    filt3DMask(:,imSlice,:)=imfill(squeeze(filt3DMask(:,imSlice,:)),'holes');
end

end

 filt3DMask=imerode(filt3DMask,se);
 filt3DMask=imerode(filt3DMask,true(1,1,5));
 filt3DMask(pedMask)=0;

for imSlice=size(filt3DMask,3)-2:size(filt3DMask,3);
    filt3DMask(:,:,imSlice)=0;
end% filt3DMask=imopen(filt3DMask,se);

% hits=(sum(sum(filt3DMask,1),2))>500;
% filt3DMask=bsxfun(@times,filt3DMask,hits);
%%

% filt3DMask= imdilate(filt3DMask,se);
% filt3DMask=imfill(filt3DMask,'holes');
% filt3DMask=imerode(filt3DMask,true(9,9,9));


%F=imerode(filt3DMask,ones(1,1,5));

% if number of pixels in mask is too small, replace with sphere centered at
% middle.
if 0
%     if mean(filt3DMask(:))>0.005; 
    CLimdist=bwdist(filt3DMask);
    F=isosurface(CLimdist,0);
    F.vertices=F.vertices(:,[2,1,3]);
    
    %F=triinterp2(F);
    
    
    %%
    %cline_para.upFactor=1;
    %F=triInterp2(F);
    % if size(F.vertices,1)<1000
    % F=surfaceResample(F,4,2);
    % F=triInterp2(F);
    % %F=surfaceResample(F,1.6,.8);
    % end
    F=surfaceResample(F,cline_para.upperLimit,cline_para.lowerLimit);
    
    F=faceNormals(F); %sort Face ordering so that normals will all point the same way.
else
    centerPix=round(image_para.imsize/2);
    filt3DMask=zeros(size(hyper_stack));
    filt3DMask(centerPix(1),centerPix(2),centerPix(3))=1;
    CLimdist=bwdist(filt3DMask);
    F=isosurface(CLimdist,min(centerPix/3));
    F.vertices=F.vertices(:,[2,1,3]);
    F=faceNormals(F); %sort Face ordering so that normals will all point the same way.
    
end

%% Notes on initial configuration
% sometimes, this initial sphere has very fine faces, maybe starting with some decimated version would be better
%% SurfaceFit
CLscale=2;

% %%% here's a place to making a starting F

%%


switch flag.isStartSphere
       case 0
   hyper2 = hyper_stack>quantile(hyper_stack(:),0.9);
   hyper2 = imdilate(hyper2,ones(5,5,5));
   hyper2 = imerode(hyper2,ones(7,7,7));
   F2 = isosurface(hyper2,0.5);
   %% subsample these surface vertices and build an alpha shape 
   divisionFactor = 5;
   randomVertices = randperm(length(F2.vertices),round(length(F2.vertices)/divisionFactor));
   
   F3 = alphaShape(F2.vertices(randomVertices,2),...
       F2.vertices(randomVertices,1),...
       F2.vertices(randomVertices,3),7);
   % Sometimes there's another object in the field of view, pick the largest object before continuing
   testThresh = 0.01;
   while F3.numRegions>1
       testThresh = testThresh*1.3;
       F3.RegionThreshold = testThresh;
   end
   %
   F4.faces = F3.boundaryFacets;
   F4.vertices = F2.vertices(randomVertices,[2,1,3]);
   F4 = faceNormals(F4);
  
%%

% ticA = tic;
   [~,vertAreaI] = triArea(F4);
%    toc(ticA);
   upperLimit = cline_para.upperLimit; lowerLimit = cline_para.lowerLimit;
   
   [F42,vertAreaI]=surfaceResample4(F4,1.5*upperLimit,0.8*lowerLimit,vertAreaI);
%    toc(ticA)
   
   
   
   F42 = faceNormals(F42);
%    %% decimate again
%       randomVertices = randperm(length(F42.vertices),...
%           round(length(F42.vertices)/2));
% 
%    F5 = alphaShape(F42.vertices(randomVertices,2),...
%        F42.vertices(randomVertices,1),...
%        F42.vertices(randomVertices,3),25);
%    % Sometimes there's another object in the field of view, pick the largest object before continuing
%    testThresh = 0.01;
%    while F5.numRegions>1
%        testThresh = testThresh*1.3;
%        F5.RegionThreshold = testThresh;
%    end
%    %
%    F6.faces = F5.boundaryFacets;
%    F6.vertices = F42.vertices(randomVertices,[2,1,3]);
%    F6 = faceNormals(F6);
   %%
   
   
   [F42,vertAreaI]=surfaceResample4(F42,1.3*upperLimit,0.9*lowerLimit,vertAreaI);
%    toc(ticA)
   F42 = faceNormals(F42);
   [F42,vertAreaI]=surfaceResample4(F42,1.2*upperLimit,0.9*lowerLimit,vertAreaI);
%    toc(ticA)
   F42 = faceNormals(F42);
   [F42,vertAreaI]=surfaceResample4(F42,1.1*upperLimit,0.9*lowerLimit,vertAreaI);
%    toc(ticA)
   %%
   F42 = faceNormals(F42);
   [F42,vertAreaI]=surfaceResample4(F42,1.0*upperLimit,0.9*lowerLimit,vertAreaI);
%    toc(ticA)
   F42 = faceNormals(F42);
   [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);

    disp(length(vertAreaI));
%     toc(ticA);
  %%
%    ticA = tic;
% 
%  upperLimit = 5; lowerLimit = 0.2
%    [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);
% disp(length(vertAreaI));
% toc(ticA);
% 
%    %%
%    ticA = tic;
% 
%  upperLimit = 5; lowerLimit = 0.3
%    [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);
% disp(length(vertAreaI));
% toc(ticA);
% 
% %%
% ticA = tic;
% 
%  upperLimit = 5; lowerLimit = 0.5
%    [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);
% disp(length(vertAreaI));
% toc(ticA);
% 
% %%
% ticA = tic;
% 
%  upperLimit = 5; lowerLimit = 0.8
%    [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);
% disp(length(vertAreaI));
% toc(ticA);
% 
% %%
% ticA = tic;
% 
%  upperLimit = 5; lowerLimit = 1.2;
%    [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);
% disp(length(vertAreaI));
% toc(ticA);
% 
% %%
% 
%  upperLimit = 5; lowerLimit = 1.5;
%    [F42,vertAreaI]=surfaceResample4(F42,upperLimit,lowerLimit,vertAreaI);
% disp(length(vertAreaI));

   %%
% clf; patch(F6,'faceColor','w','edgecolor','k');axis equal;
   %%
%    keyboard;
F = F42;
    case 2
        hyper_stack2 = cat(3,zeros(size(hyper_stack,1),size(hyper_stack,2),2),...
            hyper_stack,...
            zeros(size(hyper_stack,1),size(hyper_stack,2),2));
        F4 = isosurface(hyper_stack2,0.4);
        F4.vertices(:,3) = (F4.vertices(:,3)-2)*30/51+8;
%
   F4 = reducepatch(F4,300);
   F5 = alphaShape(F4.vertices(:,1),F4.vertices(:,2),F4.vertices(:,3),13);
   F6.faces = F5.boundaryFacets;
   F6.vertices = F5.Points(:,[2,1,3]);

%    clf
%     patch(F6,'FaceColor','w','edgecolor','k'); axis equal
    F = faceNormals(F6);
    case 3
        %%
        hyper_stack2 = cat(3,zeros(size(hyper_stack,1),size(hyper_stack,2),2),...
            hyper_stack,...
            zeros(size(hyper_stack,1),size(hyper_stack,2),2));
        
%         hyper_stack3 = hyper_stack2>0.05;
        hyper_stack3 = imerode(hyper_stack2,strel('ball',5,4));
        hyper_stack3 = normalizeRange(hyper_stack3);
        
        F4 = isosurface(hyper_stack3,0.4);
%         F4.vertices(:,3) = (F4.vertices(:,3)-2)*30/51+8;
%

    % reduce to half the number of faces
   F4 = reducepatch(F4,0.5);
   % switch to an alpha hull
   F5 = alphaShape(F4.vertices(:,1),F4.vertices(:,2),F4.vertices(:,3),5);
   F6.faces = F5.boundaryFacets;
   F6.vertices = F5.Points(:,[2,1,3]);
    
   %%
   clf
   patch(F6,'FaceColor','w','edgecolor','k'); axis equal
   
   %%
%    dt1 = delaunayTriangulation(F6.vertices);
%    tetramesh(dt1,'FaceColor','blue','FaceAlpha',0.3);

   %%
%%   
F = faceNormals(F6);
%  [~,vertAreaI] = triArea(F6);
%    toc(ticA);
%    upperLimit = cline_para.upperLimit; lowerLimit = cline_para.lowerLimit;
%    [F62,vertAreaI]=surfaceResample4(F4,15,0.5,vertAreaI);    
   %%
%    F = faceNormals(F62);
case 4
        hyper_stack2 = cat(3,zeros(size(hyper_stack,1),size(hyper_stack,2),2),...
            hyper_stack,...
            zeros(size(hyper_stack,1),size(hyper_stack,2),2));
        F4 = isosurface(hyper_stack2,0.4);
        % compress the z Height a little in the starting shape
        zMean = mean(F4.vertices(:,3));
        F4.vertices(:,3) = ((F4.vertices(:,3))-zMean)*0.75+zMean;
%
   F4 = reducepatch(F4,0.3);
   
   F5 = alphaShape(F4.vertices(:,1),F4.vertices(:,2),F4.vertices(:,3),8);
   F6.faces = F5.boundaryFacets;
   F6.vertices = F5.Points(:,[2,1,3]);

%    clf
%     patch(F6,'FaceColor','w','edgecolor','k'); axis equal
    F = faceNormals(F6);
    
    case 6
        %%
        hyper_stack2 = cat(3,zeros(size(hyper_stack,1),size(hyper_stack,2),15),...
            hyper_stack,...
            zeros(size(hyper_stack,1),size(hyper_stack,2),15));
%        hyper_stack2 = hyper_stack;
kern2 = ones(15,15,15);
hyper_stack3 = (imclearborder(hyper_stack2>0.3));
hyper_stack3 = convn(hyper_stack3,kern2./sum(kern2(:)),'same');
hyper_stack3 = double(hyper_stack3>0.05);

hyper_stack2 = convn(hyper_stack2,kern2./sum(kern2(:)),'same');
% now we have a degraded version of the original image

% 
%        hyper_stack2 = bpass3_jn(hyper_stack2,[3],[35]);
%        hyper_stack2(:,:,1:14) = [];
%        hyper_stack2 = hyper_stack2>0.09;
       
%        hyper_stack2 = imfill(hyper_stack2,'holes');
%               SliceBrowser(hyper_stack2);

        %%
        F4 = isosurface(hyper_stack2.*hyper_stack3,0.15);
        % compress the z Height a little in the starting shape
        zMean = mean(F4.vertices(:,3));
        F4.vertices(:,3) = ((F4.vertices(:,3))-15);%zMean)*0.75+zMean;
%
%    F4 = reducepatch(F4,0.8);
   
   F5 = alphaShape(F4.vertices(:,1),F4.vertices(:,2),F4.vertices(:,3),8);
   F6.faces = F5.boundaryFacets;
   F6.vertices = F5.Points(:,[2,1,3]);
% 
%    clf
%     patch(F6,'FaceColor','w','edgecolor','k'); axis equal
    F = faceNormals(F6);
    case 7
         %%
        hyper_stack2 = cat(3,zeros(size(hyper_stack,1),size(hyper_stack,2),15),...
            hyper_stack,...
            zeros(size(hyper_stack,1),size(hyper_stack,2),15));
%        hyper_stack2 = hyper_stack;
kern2 = ones(5,5,5);
hyper_stack3 = (imclearborder(hyper_stack2>0.3));
hyper_stack3 = convn(hyper_stack3,kern2./sum(kern2(:)),'same');
hyper_stack3 = double(hyper_stack3>0.05);

hyper_stack2 = convn(hyper_stack2,kern2./sum(kern2(:)),'same');
% now we have a degraded version of the original image

% 
%        hyper_stack2 = bpass3_jn(hyper_stack2,[3],[35]);
%        hyper_stack2(:,:,1:14) = [];
%        hyper_stack2 = hyper_stack2>0.09;
       
%        hyper_stack2 = imfill(hyper_stack2,'holes');
%               SliceBrowser(hyper_stack2);

        %
        F4 = isosurface(hyper_stack2,0.4);
%         F4 = isosurface(hyper_stack3,0.8);
%         % compress the z Height a little in the starting shape
%         zMean = mean(F4.vertices(:,3));
%         F4.vertices(:,3) = ((F4.vertices(:,3))-15);%zMean)*0.75+zMean;
% %
%    F4 = reducepatch(F4,0.8);
   
   F5 = alphaShape(F4.vertices(:,1),F4.vertices(:,2),F4.vertices(:,3),20);
   F6.faces = F5.boundaryFacets;
   F6.vertices = F5.Points(:,[2,1,3]);
    
%    F6 = reducepatch(F4,0.9);
     
   [~,vertAreaI] = triArea(F6);
%       F6 = reducepatch(F6,(round(sum(vertAreaI)))/5);

%    toc(ticA);
   upperLimit = cline_para.upperLimit; lowerLimit = cline_para.lowerLimit;
   %

   upperLimit = 2;
   lowerLimit = 0.5;
   currentUpper = 7*upperLimit;
   currentLower = lowerLimit/7;

   counter1= 1;
   F62 = F6;
   ticB = tic;
   while currentUpper>upperLimit
       currentUpper = currentUpper/1.5;
       currentLower = currentLower*1.5;
       [F62,vertAreaI]=surfaceResample4(F62,currentUpper,currentLower,vertAreaI);
       F62 = faceNormals(F62);

   counter1 = counter1+1;
   disp(counter1);
   end
   toc(ticB)
   
   
   %    toc(ticA)
   %
   
   
   F62 = faceNormals(F62);
   F6 = F62;
   F6.vertices = F6.vertices(:,[2,1,3]);
   F6.vertices(:,3)= F6.vertices(:,3) - 15; % number of zeros planes added earlier
   %
   clf
    patch(F6,'FaceColor','w','edgecolor','k'); axis equal
    %
    F = faceNormals(F6);
    case 8
        %%
        hyper_stack2 = cat(3,zeros(size(hyper_stack,1),size(hyper_stack,2),15),...
            hyper_stack,...
            zeros(size(hyper_stack,1),size(hyper_stack,2),15));
%        hyper_stack2 = hyper_stack;
kern2 = ones(15,15,15);
hyper_stack3 = (imclearborder(hyper_stack2>0.3));
hyper_stack3 = convn(hyper_stack3,kern2./sum(kern2(:)),'same');
hyper_stack3 = double(hyper_stack3>0.05);

hyper_stack2 = convn(hyper_stack2,kern2./sum(kern2(:)),'same');
% now we have a degraded version of the original image

% 
%        hyper_stack2 = bpass3_jn(hyper_stack2,[3],[35]);
%        hyper_stack2(:,:,1:14) = [];
%        hyper_stack2 = hyper_stack2>0.09;
       
%        hyper_stack2 = imfill(hyper_stack2,'holes');
%               SliceBrowser(hyper_stack2);

        %%
        F4 = isosurface(hyper_stack2.*hyper_stack3,0.15);
        % compress the z Height a little in the starting shape
        zMean = mean(F4.vertices(:,3));
        F4.vertices(:,3) = ((F4.vertices(:,3))-15);%zMean)*0.75+zMean;
%
%    F4 = reducepatch(F4,0.8);
   
   F5 = alphaShape(F4.vertices(:,1),F4.vertices(:,2),F4.vertices(:,3),18);
   F6.faces = F5.boundaryFacets;
   F6.vertices = F5.Points(:,[2,1,3]);
% 
%    clf
%     patch(F6,'FaceColor','w','edgecolor','k'); axis equal
    F = faceNormals(F6);
    
    otherwise
        
end
%%
[F,eg]=activeSurfaceFitTri3(F,hyper_stack,cline_para,image_para,flag);

%%
F2=triInterp2(F);
F2=triInterp2(F2);
membrane2=coord2image3d(F2,image_para.imsize,2,1);
membrane2=imclose(membrane2,true(5,5,5));
membrane2=imdilate(membrane2,true(5,5,5));

mfill=imfill(double(membrane2~=0),'holes');
mfill=imerode(mfill,true(5,5,5));
membrane2=coord2image3d(F2,image_para.imsize,CLscale,1);
%%

%%
Fcurvature=tricurv_v02(F.faces,F.vertices);

if any(isnan(Fcurvature.k1));
    warning('CellShapeDetectorTri:errorInCurvature','nan curvature detected');
    
end


Fcurvature.k1=triSmooth(F.faces,Fcurvature.k1/image_para.nm_per_pixel);
Fcurvature.k2=triSmooth(F.faces,Fcurvature.k2/image_para.nm_per_pixel);
Fcurvature.km=triSmooth(F.faces,Fcurvature.km/image_para.nm_per_pixel);
Fcurvature.kg=triSmooth(F.faces,Fcurvature.kg/image_para.nm_per_pixel^2);


F.vertices=F.vertices*image_para.nm_per_pixel;

Fshape.Area=triArea(F);
%Fshape.Length=max(s)*image_para.nm_per_pixel;
Fshape.Volume=triVolume(F);



%% centerline
Vmat=vertDistMat(F2);
gdist=graphallshortestpaths(Vmat);
[pnt1,pnt2]=find(gdist==max(gdist(:)),1,'first');


%%
% centerline updates to find poles

Vmat=vertDistMat(F);
% working on the surface, find all shortest paths
gdist=graphallshortestpaths(Vmat);

% start estimating endpoints as max geodesic
[pnt1,pnt2]=find(gdist==max(gdist(:)),1,'first');

% look near each of these and build a neighbourhood/patch that increases in gaussian curvature
[F2Vmat,invTri]=Face2Vert(F.faces);

noiseLvl = 0.15*(max(Fcurvature.kg(:))-quantile(Fcurvature.kg(:),.2));

%% start at pnt1 
isAnotherRound = true;
% create 'sparse' array to define surface region
currentPatch = sparse(size(F2Vmat,2),1);
% start at a point
currentPatch(pnt1) = true;
% vertex to vertex connectivity
Vert2Vert = F2Vmat'*F2Vmat;
Vert2Vert(Vert2Vert>0) = 1;
% only keep gaussian curvatures that are higher (with some noise allowed)
threshGCurv = Fcurvature.kg>(Fcurvature.kg(pnt1)-noiseLvl);
patchSize = 1; % counts vertices
while isAnotherRound
    % increase size of patch to include one more set of neighbors
    currentPatch = logical((currentPatch'*Vert2Vert).*threshGCurv')';
    % check to see if there were any new neighbors added
    isAnotherRound = nnz(currentPatch)-patchSize;
    patchSize = nnz(currentPatch);
end

% look up which faces are included in this new surface patch
facesKept = invTri(currentPatch,:);
facesKept = unique(facesKept);
facesKept(facesKept==0) = [];
% build a new patch style object with these
F3 = F;
F3.faces = F3.faces(facesKept,:);


pnt1b = pnt1;
[maxCurv1,pnt1] = max(currentPatch.*Fcurvature.kg);


% 
% % 
% % %find all edges of the current patch, boundary edges appear only once. 
% % E=[F3.faces(:,[1,3]);F3.faces(:,2:3);F3.faces(:,1:2)];
% % E=sort(E,2);
% % [E,ia,ib]=unique(E,'rows');
% % boundaryEdge=accumarray(ib,ones(1,length(ib))')==1;
% % patchFreeEdgeVert=E(boundaryEdge,:);
% % patchFreeEdgeVert=unique(patchFreeEdgeVert(:));
% % patchAllVert = unique(F3.faces(:));
% 
% 
% 
% 
% 
% 
% 
% % % build a new triangulation of just the free boundary or "edge" of the
% % % patch
% % TR = triangulation(F3.faces,F3.vertices);
% % [FF.faces, FF.vertices] = freeBoundary(TR);
% % % convert back into numbering scheme of F3
% % [~,F3label] = ismember(FF.vertices,F3.vertices,'rows');
% % 
% % FF.vertices = F3.vertices;
% % FF.faces = F3label(FF.faces);
% % 
% % % pull out all the vertices of interest
% % patchFreeEdgeVert = unique(FF.faces(:));
% % patchAllVert = unique(F3.faces(:));
% 
% 
% 
% % patchDistances = gdist(patchAllVert,patchFreeEdgeVert);
% % totalDistance = sum(patchDistances,2);
% % [~,maxId] = min(var(gdist(patchAllVert,patchFreeEdgeVert),1,2));
% % pnt1b = pnt1;
% % pnt1 = patchAllVert(maxId);
% 
%

%
% display this resulting patch and pole
if cline_para.show_flag
figure;
hold on;
patch(F,'facevertexCdata',Fcurvature.kg,'facecolor','interp','facealpha',0.2,'edgecolor','none');
patch(F3,'facevertexCdata',Fcurvature.kg,'facecolor','interp','facealpha',0.6,'edgecolor','black','edgealpha',0.3);
scatter3(F.vertices(:,1),F.vertices(:,2),F.vertices(:,3),4,'kx');
%patch(FF,'facevertexCdata',Fcurvature.kg,'facecolor','none','facealpha',0.0,'edgecolor','black');
scatter3(F.vertices(pnt1b,1),F.vertices(pnt1b,2),F.vertices(pnt1b,3),250,'w+');
scatter3(F.vertices(pnt1,1),F.vertices(pnt1,2),F.vertices(pnt1,3),250,'k*');

axis equal;
end

%% start at pnt2
isAnotherRound = true;
% create 'sparse' array to define surface region
currentPatch = sparse(size(F2Vmat,2),1);
% start at a point
currentPatch(pnt2) = true;
% vertex to vertex connectivity
Vert2Vert = F2Vmat'*F2Vmat;
Vert2Vert(Vert2Vert>0) = 1;
% only keep gaussian curvatures that are higher (with some noise allowed)
threshGCurv = Fcurvature.kg>(Fcurvature.kg(pnt2)-noiseLvl);
patchSize = 1; % counts vertices
while isAnotherRound
    % increase size of patch to include one more set of neighbors
    currentPatch = logical((currentPatch'*Vert2Vert).*threshGCurv')';
    % check to see if there were any new neighbors added
    isAnotherRound = nnz(currentPatch)-patchSize;
    patchSize = nnz(currentPatch);
end

% look up which faces are included in this new surface patch
facesKept = invTri(currentPatch,:);
facesKept = unique(facesKept);
facesKept(facesKept==0) = [];
% build a new patch style object with these
F3 = F;
F3.faces = F3.faces(facesKept,:);

pnt2b = pnt2;
[maxCurv2,pnt2] = max(currentPatch.*Fcurvature.kg);



% 
% %find all edges of the current patch, boundary edges appear only once. 
% E=[F3.faces(:,[1,3]);F3.faces(:,2:3);F3.faces(:,1:2)];
% E=sort(E,2);
% [E,ia,ib]=unique(E,'rows');
% boundaryEdge=accumarray(ib,ones(1,length(ib))')==1;
% patchFreeEdgeVert=E(boundaryEdge,:);
% patchFreeEdgeVert=unique(patchFreeEdgeVert(:));
% patchAllVert = unique(F3.faces(:));
% 
% 
% 
% % % build a new triangulation of just the free boundary or "edge" of the
% % % patch
% % TR = triangulation(F3.faces,F3.vertices);
% % [FF.faces, FF.vertices] = freeBoundary(TR);
% % % convert back into numbering scheme of F3
% % [~,F3label] = ismember(FF.vertices,F3.vertices,'rows');
% % 
% % FF.vertices = F3.vertices;
% % FF.faces = F3label(FF.faces);
% % 
% % % pull out all the vertices of interest
% % patchFreeEdgeVert = unique(FF.faces(:));
% % patchAllVert = unique(F3.faces(:));
% 
% 
% 
% % patchDistances = gdist(patchAllVert,patchFreeEdgeVert);
% % totalDistance = sum(patchDistances,2);
% [~,maxId] = min(var(gdist(patchAllVert,patchFreeEdgeVert),1,2));
% pnt2b = pnt2;
% pnt2 = patchAllVert(maxId);
%
if cline_para.show_flag 
patch(F3,'facevertexCdata',Fcurvature.kg,'facecolor','interp','facealpha',0.6,'edgecolor','black','edgealpha',0.3);
scatter3(F.vertices(:,1),F.vertices(:,2),F.vertices(:,3),4,'kx');
%patch(FF,'facevertexCdata',Fcurvature.kg,'facecolor','none','facealpha',0.0,'edgecolor','black');
scatter3(F.vertices(pnt2b,1),F.vertices(pnt2b,2),F.vertices(pnt2b,3),250,'w+');
scatter3(F.vertices(pnt2,1),F.vertices(pnt2,2),F.vertices(pnt2,3),250,'k*');
 
end



%% build centerline along geodesic between these "poles"

% check to make sure that the two endpoints are different
if isequal(pnt1,pnt2);
    warning('CellShapeDetector:TriFolder:IdenticalPoles',...
        'Both poles are found at the same point. Reverting to initial pole guess');
    pnt1 = pnt1b;
    pnt2 = pnt2b;
end

% find the geodesic from pnt1 to pnt2
[dist, path, ~]=graphshortestpath(Vmat,pnt1,pnt2);
% which points are included on the geodesic?
pathpts=zeros(1,length(F.vertices));
pathpts(path)=1;
% find the geodesic to each of the end poles
d1=gdist(pnt1,:);
d2=gdist(pnt2,:);
% calculate the fraction of the distance travelled
L=d1./(d1+d2);
% divide into 10 "hoops" based on the fraction of the way along the
% centerline these vertices are
L2=round(L*10);

% find the center of mass of each of these hoops
CL=[accumarray(L2'+1,F.vertices(:,1),[],@mean),...
    accumarray(L2'+1,F.vertices(:,2),[],@mean),...
    accumarray(L2'+1,F.vertices(:,3),[],@mean)];

% this curve's contour variable is called s
s=[0,cumsum(sqrt(diff(CL(:,1)).^2+diff(CL(:,2)).^2+diff(CL(:,3)).^2))']';
 
    CL=[interp1(s,CL(:,1),0:image_para.nm_per_pixel:max(s),'spline')',...
        interp1(s,CL(:,2),0:image_para.nm_per_pixel:max(s),'spline')',...
        interp1(s,CL(:,3),0:image_para.nm_per_pixel:max(s),'spline')'];

CL(1,:)=F.vertices(pnt1,:);
CL(end,:)=F.vertices(pnt2,:);


%use membrane to create distance energy that the Centerline will follow,
%keeping in mind to take account for scaling of the membraneimage. 

        M3=bwdist(membrane2).^2;
        
        CLind=round(CL*CLscale/image_para.nm_per_pixel);
        CLind=sub2ind(size(M3),CLind(:,1),CLind(:,2),CLind(:,3));
        M3=M3/max(M3(CLind));
        CL = ActiveContourFit(M3, cline_para, CL*CLscale/image_para.nm_per_pixel);
        CL=CL/CLscale*image_para.nm_per_pixel;
        
        s=[0,cumsum(sqrt(diff(CL(:,1)).^2+diff(CL(:,2)).^2+diff(CL(:,3)).^2))']';
        si=0:.1:max(s);
        %     CL=[interp1(s,CL(:,1),si,'spline')',...
        %         interp1(s,CL(:,2),si,'spline')',...
        %         interp1(s,CL(:,3),si,'spline')'];
        %
        % CL(1,:)=F.vertices(pnt1,:);
        % CL(end,:)=F.vertices(pnt2,:);
        
        %% Calculate Radii, and nearest Radii point
        
        
        %%
        % theta=zeros(length(F2.vertices),1);
        % for i=2:length(path)-1
        % Ip1=I(path(i-1));
        % Ip2=I(path(i+1));
        %     subpts=(I>Ip1&I<Ip2)'; %vertices at some Lrange
        %  %subpts=subpts*subpts';
        %  Vmattemp=Vmat2;
        %     Vmattemp(~subpts,~subpts)=0;
        %     subdist=graphallshortestpaths(Vmattemp);
        %     rloop=subdist(:,path);
        %     rloop(rloop==0)=Inf;
        %     rloop=min(rloop,[],2);
        %     rloop=rloop/sax(rloop(~isinf(rloop)));
        %     theta(~isinf(rloop))=rloop(~isinf(rloop));
        % end
        %
        % Fcoord.theta=theta;
        % Fcoord.length=I;
        
        
      
        %% calculate tangent, normal, binormal of the centerline

        [Tv,Bv,Nv]=tbnVector(CL);
        Bv=[Bv(1,:);Bv];
        Nv=[Nv(1,:);Nv];
        
        % calculate distance to the centerline

        [R,I2] = pdist2(CL,F.vertices,'euclidean','Smallest',1);
        % projection along the centerline that these vertices are closest
        % to
        I=s(I2);
        % two unused variables??  
        F2Vmat=Face2Vert(F.faces);
        Vmat2=vertDistMat(F,CL,'circular');
        % moving into cylindrical coordinates by drawing vectors from the
        % centerline to the surface
        Fcyl=F.vertices-CL(I2,:);
        % project these vectors onto the normal and binormal
        xcyl=dot(Fcyl,Nv(I2,:),2);
        ycyl=dot(Fcyl,Bv(I2,:),2);
        % use this as a measure of the unwrap angle theta
        theta=atan2(ycyl,xcyl);
        Fcoord.theta=theta;
        Fcoord.L=I';%;*image_para.nm_per_pixel;
        Fcoord.r=sqrt(xcyl.^2+ycyl.^2);
        Fcoord.CL=CL;%*image_para.nm_per_pixel;
        
        
        
        %         CL=[mean(x)' mean(y)' mean(z)'];
        %         M3=bwdist(~mfill).^2;
        %         M3=M3/max(M3(:))*CLscale;
        %         CL = ActiveContourFit(M3, cline_para, xyzs*CLscale);
        
        
        %%
        
        
        %% Calculate surface area, volume and centerline curvature;
% Volume=surfVolume(xf,yf,zf);
% K=Curvature(CL,10);
% T=Twist(CL,10);

%% Create Fluorescent map by interpolating

if flag.F1==1
    
    
    for tt=1:image_para.Fstack_t_size
        % if theres a Fluorescent stack, open it and scale it so that the
        %axis are equal and have its center be the same as the membrane stack
        F1=zeros(row,col,image_para.Fstack_z_size);
        for k = 1:image_para.Fstack_z_size
            if flag.reverse==0
                                      planeIndex =  image_para.stack_z_size*image_para.stack_t_size+(tt-1)*image_para.Fstack_z_size+k;
                        
                    else
                        planeIndex = (tt-1)*image_para.Fstack_z_size+k;
                        
                    end
                    F1(:,:,k) = double(imread([p,'.tif'], 'tif',planeIndex))+...
                        intensity_correction;
                end
                %%% ----- filter shift patch
                % check for filter information
                if isMetaData
                    % background subtract, shiftImgFilter, background add
                    backgndTmp = quantile(F1(:),0.005);
                    F1 = F1-backgndTmp;
                    for k2 = 1:image_para.Fstack_z_size
                        if flag.reverse==0
                            planeIndex =  image_para.stack_z_size*image_para.stack_t_size+(tt-1)*image_para.Fstack_z_size+k2;
                            
                        else
                            planeIndex = (tt-1)*image_para.Fstack_z_size+k2;
                            
                        end
                        
                        filterString = totalImageInfo{planeIndex}.filterID;
                        filterIdx = str2double(filterString(end-1));
                        if isnan(filterIdx)
                            filterIdx = 3; % no shifting of image
                        end
                        
                        F1(:,:,k2) = shiftImgFilter(F1(:,:,k2),filterIdx);
                        
                    end
                    F1 = F1+backgndTmp;
                end
                %%% ----- end filter shift patch        
        F1=image_resize(F1,size(F1,1),...
    size(F1,2),round(image_para.ZF_factor*size(F1,3)));
        
        % This will show the fit outlines on the fluorescent channel
       
        %if the is suspected motion between to 2 channels in addition to
        %the correction from the filters, find the correlation and shift
        %the image;
        
        if isfield(flag,'F1fix')
            if flag.F1fix == 1
                F1=F1-median(F1(pedMask));
                F1(F1<0)=0;
                AAA=convnfft(hyper_stack.*~pedMask,...
                    F1(end:-1:1,end:-1:1,end:-1:1).*~pedMask);
                [xfix,yfix,zfix]=ind2sub(size(AAA),find(AAA==max(AAA(:))));
                xfix=xfix-ceil(size(AAA,1)/2);
                yfix=yfix-ceil(size(AAA,2)/2);
                zfix=zfix-ceil(size(AAA,3)/2);
                
                
                [Xmap,Ymap,Zmap]=ndgrid(1:size(F1,1),1:size(F1,2),1:size(F1,3));
                Xmap=Xmap-xfix;
                Ymap=Ymap-yfix;
                Zmap=Zmap-zfix;
                F1=interp3(F1,Ymap,Xmap,Zmap,'nearest',0);
            elseif flag.F1fix == 2
                F1 = bpass3_jn(F1,[2],[11,11,11]);
            end
            
        end
        
        %This is set up to measure the intensities at the membrane for each
        %point.
        Ff=interp3(F1,F.vertices(:,2)/image_para.nm_per_pixel...
            ,F.vertices(:,1)/image_para.nm_per_pixel...
            ,F.vertices(:,3)/image_para.nm_per_pixel,'*linear');
        
        %max radial projection
        
        r_steps=0.5:.05:2;  %radius range
        
        r_steps=reshape(r_steps,1,1,[]);
        

                r_vectors=(F.vertices/image_para.nm_per_pixel-CL(I2,:));
                r_vectors=bsxfun(@times,r_vectors,r_steps);
                r_vectors=bsxfun(@plus,r_vectors,CL(I2,:));
                r_vectors=round(r_vectors);
                r_vectors(r_vectors<=0)=1;
                r_vectors=bsxfun(@min,r_vectors,size(F1));
                
%F2=interp3(F1,r_vectors(:,1,:),r_vectors(:,2,:),r_vectors(:,3,:),'*linear');
F2=F1(sub2ind(size(F1),r_vectors(:,1,:),r_vectors(:,2,:),r_vectors(:,3,:)));
F2=squeeze(F2);

                Ffr=max(F2,[],2);
                %Ftest(k,i,:)=F2;
                intensityMap.patch_at_surface=Ff;
                intensityMap.patch_Rprojected=Ffr;
                %% Ben thinks this must be a typo, I is the centerline coordinate
                intensityMap.FM_at_surface=I;
                
                
            end
        end
        
        
        %% Save files
        if flag.F1==1
            if image_para.Fstack_t_size>1
                FF=Ff;
                for tt=1:image_para.Fstack_t_size
                    p2=[p num2str(tt,'%02d') ptail ];
                    Ff=FF(:,:,tt);
                    
                    save(p2,'p2','F','Fcurvature','Fshape')
                    
                end
                
            else
                p2=[p,ptail];
            
            
            save(p2,'p2','F','Fcurvature','Fshape','eg','Ff','Ffr','Fcoord','image_para','cline_para')
            
            end
        end
        
        if flag.F1==0
            p2=[p,ptail];
            
            save(p2,'p2','F','Fcurvature','Fshape','eg','Fcoord')
            
        end
        
        
        
        %       clear xpe ype zpe xpef ypef zpef xp2ef yp2ef zp2ef xps yps zps xp2sf yp2sf zp2sf xpsf ypsf zpsf
        %       clear xp2 yp2 zp2 xp2e yp2e zp2e xp2s yp2s zp2s Ftest cl centerline
        
        %%
        if Fflag2==1 %%this is for DNA stain or other signal on interior of cell
            warning('cellShape:untestedParameter:Fflag2','Fflag2 set to true. This is an untested mode.');
            
            F=zeros(row,col,image_para.Fstack_z_size);
            
            for k = 1:Fstack_z_size
                F(:,:,k) = double(imread([p,'.tif'], 'tif',image_para.stack_z_size*image_para.stack_t_size*2+k)+intensity_correction);
            end
            F=image_resize(F,size(F,1),size(F,2),round(image_para.ZF_factor*size(F,3)));
            [X,Y,Z]=meshgrid(1:col,1:row,1:image_para.stack_z_size);
            [FX,FY,FZ]=meshgrid(1:col,1:row,1:size(F,3));
            F1=interp3(FX,FY,FZ-size(F,3)/2,F,X,Y,Z-image_para.stack_z_size/2,'nearest');
            
            for k=1:image_para.stack_z_size %repeat the ends of the fluorescent stack so all points beyond the end will look like the end. This is for extrapolation
                if k<image_para.stack_z_size/2 && all(all(isnan(F1(:,:,k))))
                    F1(:,:,k)=F(:,:,1);
                end
                if k>image_para.stack_z_size/2 && all(all(isnan(F1(:,:,k))))
                    F1(:,:,k)=F(:,:,end);
                end
            end
            clear F
            %This is set up to measure the average F intensity inside the
            %contour of each slice
            
            Ff2=zeros(size(xys,3),1);
            for iSlice=1:size(xys,3)
                F=interp3(X,Y,Z,F1,yslice(:,:,iSlice),xslice(:,:,iSlice),zslice(:,:,iSlice),'linear',NaN);
                bwmask=poly2mask(xys(:,1,iSlice),xys(:,2,iSlice),2*window+1,2*window+1);
                bwmask(isnan(F))=0;
                Ff2(iSlice)=mean(F(bwmask));
                %             imagesc(F(:,:,i));axis equal
                %             hold on
                %             plot(xys(:,1,i),xys(:,2,i))
                %             pause(.1)
                %
            end
            
            
            %             Ff2=interp3(X,Y,Z,F1,centerline(:,2)/nm_per_pixel,...
            %                 centerline(:,1)/nm_per_pixel,centerline(:,3)/nm_per_pixel);
            clear F FX FY FZ
            
            
        else Ff2=[];
        end
        
        disp(['Finished ',stack_name,'. Elapsed time ',num2str(toc(ticA)),' seconds.']);
    catch ME
        save(fullfile(folderName,[stack_name,'TRIerrorCaught.mat']));
%         rethrow(ME);
%         keyboard;
        disp(['Error in file ',stack_name,'. Skipping and proceeding with next.',...
            'Elapsed time ',num2str(toc(ticA)),' seconds.']);
    end
    
end
% end loop over files -----

%% Save the relevent variables


% disp('done');

%% plot
%surf(x0(1:20,90:end),y0(1:20,90:end),z0(1:20,90:end));
%
% surf(xf,yf,zf,gc)
% axis equal
% caxis([min(min(gc))/5,-min(min(gc))/5]);
