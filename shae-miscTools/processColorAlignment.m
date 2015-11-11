%% select folder to process
[FolderName] = uigetdir();
%
[fileList] = dir([FolderName,filesep,'*.tif']);

%% find bead positions in these stacks
fovID = nan(1,size(fileList,1));
fovID(1) = 1;
for iFile = 1:size(fileList,1);
%     progressbar((iFile-1)/size(fileList,1));
    %     tic;
    imgInformation = imfinfo(fullfile(FolderName,fileList(iFile).name));
    nPlanes = size(imgInformation,1);
    imgStack = nan([imgInformation(1).Width,imgInformation(1).Height,nPlanes]);
    for iiPlane = 1:nPlanes;
        imgStack(:,:,iiPlane) = imread(fullfile(FolderName,fileList(iFile).name),'Index',iiPlane,'info',imgInformation);
    end
    
    
    % peak finding parameters
    options.method = 1; % (1) centroid (2) circular symmetry
    options.interactive = 0; % don't display anything
    peaksStored = [];
    pkfndThresh = 7;
    %     pkfndThresh = 250;
    pkfndSize = 11;
    centfindSize = 11;
    bpassSize = 5;
    bpassNoise = 1;
    %     bpassSize = 0;
    %     bpassNoise = 0;
    
    % process each plane
    for iiPlane = 1:nPlanes;
        % filter for zero background
        tmpPlane = bpass(imgStack(:,:,iiPlane),bpassNoise,bpassSize);
        
        % initial peak estimate
        peaks = pkfnd(tmpPlane,pkfndThresh,pkfndSize);
        
        % refine peak estimate with center of mass
        peaks = centfind(tmpPlane,peaks,pkfndSize,bpassSize,options);
        
        % slow method of building peaks array
        peaks(:,5) = iiPlane;
        peaksStored = [peaksStored; peaks];
    end
    
    %% match into "trajectories", rather easy for the case of stationary beads
    param = [];
    param.mem = 2; % allows particles to 'blink'
    minFrames = round(nPlanes*0.8); % only interested in ones that are always found
    param.good = minFrames ;
    trackMaxDisp = 4; % upper bound estimate size of jumps between frames
    trackIt = track(peaksStored(:,[1,2,5]),trackMaxDisp,param);
    
    %% separate into trajectories
    traj = struct();
    stdPosSave = [];
    meanPosSave = [];
    
    for iiTraj = 1:trackIt(end,4)
        traj(iiTraj).positions = trackIt(trackIt(:,4)==iiTraj,1:3);
        % reject noise
        if size(traj(iiTraj).positions,1)<minFrames
            continue
        end
        traj(iiTraj).meanPos = mean(traj(iiTraj).positions);
        traj(iiTraj).stdPos = std(traj(iiTraj).positions);
        stdPosSave(iiTraj,1:2) = traj(iiTraj).stdPos(1:2);
        meanPosSave(iiTraj,1:2) = traj(iiTraj).meanPos(1:2);
        meanPosSave(iiTraj,3) = size(traj(iiTraj).positions,1);
    end
    
    % remove uninteresting entries
    traj(meanPosSave(:,3)==0) = [];
    stdPosSave(meanPosSave(:,3)==0,:) = [] ;
    meanPosSave(meanPosSave(:,3)==0,:) = [];
    
    %% display image
%     figure(gcf);
%     imshow(mean(imgStack,3),[],'InitialMagnification','fit');
%     hold on;
%     scatter(meanPosSave(:,1),meanPosSave(:,2),'rx');
%     title(fileList(iFile).name,'Interpreter','none');
%     
    %     %% ask user if they want to continue
    %     answer = questdlg('Continue with these settings?');
    %
    %     if strcmp(answer,'Yes')
    %         %% save positions
    save(fullfile(FolderName,[fileList(iFile).name,'.mat']),'traj','meanPosSave','stdPosSave');
    %     else
    %         break
    %     end
    
%     %% ask user if it is a new field
        if iFile>1    
%             answer = questdlg('Same field of view as previous?');
%     
%             if strcmp(answer,'Yes')
                fovID(iFile) = fovID(iFile-1);
%             else
%                 fovID(iFile) = fovID(iFile-1)+1;
%             end
        end
    % toc;
end
disp('Finished finding peaks');

%% load sets of positions and compare
filterID = [];
clear imageTags;
for iFile = 1:size(fileList,1);
    imageTags(iFile) = parseTags(fullfile(FolderName,fileList(iFile).name));
    filterID(iFile,1) = str2double(imageTags(iFile).filterID(end-1));
end

%% separate by field of view
% fovID = ones(size(fileList,1),1);
% fovID(2:end) = fovID(2:end)+1;
% fovID(8:end) = fovID(8:end)+1;
% fovID(20:end) = fovID(20:end)+1;
% fovID(26:end) = fovID(26:end)+1;

%% process fields of view
dtDist = 5; % keep beads if they are within this distance in two color channels
results = [];
for iiFieldOfView = 1:fovID(end);
    % for iiFieldOfView = 3;
    primaryIDX = find(and(fovID==iiFieldOfView,filterID'==3),...
        1,'first');
    if not(isempty(primaryIDX))
        load(fullfile(FolderName,[fileList(primaryIDX).name,'.mat']));
        idx2remove = [];
        %         switch primaryIDX
        %             case 4 % most are pretty bad in this field of view
        %                 idx2remove = [4];
        %             case 29
        %                 idx2remove = [5,13];
        %             case {11,22}
        %                 idx2remove = [];
        %             otherwise
        %                 idx2remove = [];
        %         end
        meanPosSave(idx2remove,:) = [];
        stdPosSave(idx2remove,:) = [];
        meanPosSaveBase = meanPosSave;
        stdPosSaveBase = stdPosSave;
        
        % some fields of view have excess particles, only use the close ones
        baseDelaunTri = DelaunayTri(meanPosSave(:,1:2));
        for iiFile = 1:size(fileList,1);
            if and(fovID(iiFile)==iiFieldOfView,iiFile~=primaryIDX)
                load(fullfile(FolderName,[fileList(iiFile).name,'.mat']));
                [idx,dist] = nearestNeighbor(baseDelaunTri,meanPosSave(:,1:2));
                meanPosSave = meanPosSave(dist<dtDist,:);
                stdPosSave = stdPosSave(dist<dtDist,:);
                % sometimes the indexing of the beads gets a little off, be
                % sure to match the beads to themselves
                meanPosSaveBaseTemp = meanPosSaveBase(idx(dist<dtDist),:);
                
                % simple shift
                diffPos = meanPosSave(:,1:2)-meanPosSaveBaseTemp(:,1:2);
                % asusme one gaussian, see if they all fit
                fitGm = gmdistribution.fit(diffPos,1);
                pdf1 = fitGm.pdf(diffPos);
                
                idx2remove = log(pdf1)<median(log(pdf1)-3);
                nnz(idx2remove)
                % after removing outliers, continue with finding mapping
                meanPosSaveBaseTemp(idx2remove,:) = [];
                meanPosSave(idx2remove,:) = [];
                
                % simple shift
                diffPos = meanPosSave(:,1:2)-meanPosSaveBaseTemp(:,1:2);
                meanDiff = mean(diffPos);  % average shift
                stdDiff = std(diffPos);
                % error associated with this mapping
                new1 = bsxfun(@minus,meanPosSave(:,1:2),meanDiff)-meanPosSaveBaseTemp(:,1:2);
                err1 = sqrt(mean(sum(new1.^2,2)));
                
                % modified affine (only scaling the same is xy)
                tform2 = cp2tform(meanPosSave(:,1:2),meanPosSaveBaseTemp(:,1:2),'affine');
                tformMat = tform2.tdata.T;
                xslide = tformMat(3,1);
                yslide = tformMat(3,2);
                scale = 0.5*(tformMat(1,1)+tformMat(2,2));
                tform2 = maketform('affine',[[scale,0,0];[0,scale,0];[xslide,yslide,1]]);
                % error associated with this mapping
                new2 = tformfwd(tform2,meanPosSave(:,1:2))-meanPosSaveBaseTemp(:,1:2);
                err2 = sqrt(mean(sum(new2.^2)));
                
                % store results
                results{iiFile,1} = num2str(filterID(iiFile));
                if err2>err1
                    results{iiFile,2} = num2str(meanDiff(1));
                    results{iiFile,3} = num2str(meanDiff(2));
                    results{iiFile,4} = '1';
                    results{iiFile,5} = num2str(err1);
                    results{iiFile,6} = num2str(err2);
                else
                    results{iiFile,2} = num2str(xslide);
                    results{iiFile,3} = num2str(yslide);
                    results{iiFile,4} = num2str(scale);
                    results{iiFile,5} = num2str(err1);
                    results{iiFile,6} = num2str(err2);
                end
                
                % include file names in results
                results{iiFile,7} = fileList(primaryIDX).name;
                results{iiFile,8} = fileList(iiFile).name;
                results{iiFile,9} = sqrt(mean(sum(stdPosSave.^2,2)));

                
%                 % display peaks found and displacements
%                 clf;
%                 scatter(meanPosSave(:,1),meanPosSave(:,2),'rx');
%                 hold on;
%                 scatter(meanPosSaveBaseTemp(:,1),meanPosSaveBaseTemp(:,2),'bo');
%                 quiver(meanPosSave(:,1),meanPosSave(:,2),diffPos(:,1),diffPos(:,2));
%                 xlim([1,512]);
%                 ylim([1,512]);
%                 axis image;
%                 title(fileList(iiFile).name,'Interpreter','none');
%                 pause();
            end
        end % loop over files with this field of view
    end % if field of view has a base image
end % loop over field of view



