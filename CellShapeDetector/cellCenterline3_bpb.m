function [xyzs,zbw,zsq]=cellCenterline3_bpb(hyper_stack,z_level,window, stack_z_size2,flag)

%% This is a new centerline detection technique
lnoise = 0.75;
lobj = 9;
threshold = 0.2;
upFactor=2;


% bandpass in 3D
filt3DStack = bpass3(hyper_stack,lnoise,lobj);
% threshold
filt3DMask = (filt3DStack./max(filt3DStack(:)))>0.2;
% find coordinates of voxels remaining
IND = find(filt3DMask);
[surfPts(:,1),surfPts(:,2),surfPts(:,3)] = ind2sub(size(filt3DMask),IND);

% clean up surface points
surfPts(surfPts(:,1)<lobj+4,:) = [];
surfPts(surfPts(:,2)<lobj+4,:) = [];
surfPts(surfPts(:,3)<lobj+4,:) = [];
surfPts(surfPts(:,1)>size(filt3DStack,1)-lobj-4,:) = [];
surfPts(surfPts(:,2)>size(filt3DStack,2)-lobj-4,:) = [];
surfPts(surfPts(:,3)>size(filt3DStack,3)-lobj-4,:) = [];

% % % display results
% % scatter3(surfPts(:,1),surfPts(:,2),surfPts(:,3));
% % axis equal;

% decimate surfPts
fracPts = 0.15;
surfPts2 = surfPts(randperm(size(surfPts,1),round(fracPts*size(surfPts,1))),:);
alphaHullRadius = 5;
[triHull, vbOutside, vbInside] = AlphaHull(surfPts2, alphaHullRadius); % works pretty well

% % % display
% % trimesh(triHull,surfPts2(:,1),surfPts2(:,2),surfPts2(:,3),'FaceColor','none');
% % axis equal;

%% now we need to convert from the triangular mesh back to cylindrical coordinates?

% convert triangulation to graph with euclidean edge distances
% calculate edge lengths
tempDist = sqrt(sum((surfPts2(triHull(:,1),:)-surfPts2(triHull(:,2),:)).^2,2));
graph12 = sparse(triHull(:,1),triHull(:,2),tempDist,length(surfPts2),length(surfPts2));
graph21 = sparse(triHull(:,2),triHull(:,1),tempDist,length(surfPts2),length(surfPts2));


tempDist = sqrt(sum((surfPts2(triHull(:,2),:)-surfPts2(triHull(:,3),:)).^2,2));
graph23 = sparse(triHull(:,2),triHull(:,3),tempDist,length(surfPts2),length(surfPts2));
graph32 = sparse(triHull(:,3),triHull(:,2),tempDist,length(surfPts2),length(surfPts2));

tempDist = sqrt(sum((surfPts2(triHull(:,3),:)-surfPts2(triHull(:,1),:)).^2,2));
graph31 = sparse(triHull(:,3),triHull(:,1),tempDist,length(surfPts2),length(surfPts2));
graph13 = sparse(triHull(:,1),triHull(:,3),tempDist,length(surfPts2),length(surfPts2));

% merge graphs together
Vmat = graph12+graph21+graph23+graph32+graph31+graph13;

% find geodesics
gdist=graphallshortestpaths(Vmat);
gdist(not(isfinite(gdist))) = nan; % not connected
[pnt1,pnt2]=find(gdist==max(gdist(:)),1,'first');
[dist, path, ~]=graphshortestpath(Vmat,pnt1,pnt2);

%% find centerline by grouping points that map onto the geodesic
% for iiPoint = 1:size(Vmat,1);
%     % pull out just the map distances between this point and the geodesic
%     tempDistance = gdist(path,iiPoint);
%     
%     % the index on the major geodesic is it's name
%     [minDist(iiPoint),pointInd(iiPoint)] = nanmin(tempDistance);
% end

% This works for all the points on the surface, those points not on the
% surface because they are interior or because they were decimated away
% should still be counted

% loop over all interior points and calculate distance to all surface
% points (shouldn't we already have this in some triangulation?)
%% group points by which main geodesic index is closest
clear pointInd tempDist
surfacePtIndx = unique(triHull(:));

for jjPoint = 1:size(Vmat,1);
    % loop over each point of the surface
    for iiPoint = 1:length(path)
        tempDist(iiPoint) = sqrt(sum((surfPts2(path(iiPoint),:)-surfPts2(jjPoint,:)).^2,2));
    end
    
    % keep the index of this one
    [~,pointInd(jjPoint)] = min(tempDist);    
end


% loop over each element of the main geodesic and average the points in
% that volume
for jjPoint = 1:length(path)
    % include a sliding window with [0.5,1,0.5] weightings
    pointsTemp1 = surfPts2(pointInd==jjPoint,:);
    pointsTemp2 = surfPts2(pointInd==jjPoint-1,:);
    pointsTemp3 = surfPts2(pointInd==jjPoint+1,:);
    CL(jjPoint,:) = mean(pointsTemp1)+0.5*mean(pointsTemp2)+0.5*mean(pointsTemp3);
    CL(jjPoint,:) = CL(jjPoint,:)/2;
    
%     pointsTemp1 = bsxfun(@minus,pointsTemp1,CL(jjPoint,:));
%     pointsTemp2 = bsxfun(@minus,pointsTemp2,CL(jjPoint,:));
%     pointsTemp3 = bsxfun(@minus,pointsTemp3,CL(jjPoint,:));
%     
    
end
% fix the end points (don't have contributions before and after)
CL(1,:) = surfPts2(path(1),:);
CL(end,:) = surfPts2(path(end),:);

% calculate an array with the individual distances along the geodesic
distGeodesic = zeros(length(path),1);
% and along the centerline contour
distCL = zeros(length(path),1);
for jjPoint = 1:length(path)-1
    distGeodesic(1+jjPoint) = gdist(path(jjPoint),path(jjPoint+1));
    distCL(1+jjPoint) = sqrt(sum((CL(jjPoint,:)-CL(jjPoint+1,:)).^2,2));
end
distGeodesic = cumsum(distGeodesic);
distCL = cumsum(distCL);

% calculate an array with distances along the centerline
xyzs(:,1)=interp1(distCL,CL(:,1),0:1:max(distCL),'spline');
xyzs(:,2)=interp1(distCL,CL(:,2),0:1:max(distCL),'spline');
xyzs(:,3)=interp1(distCL,CL(:,3),0:1:max(distCL),'spline');

% %% display
% clf;
% % trimesh(triHull,surfPts2(:,1),surfPts2(:,2),surfPts2(:,3));
% hold on;
% scatter3(surfPts2(:,1),surfPts2(:,2),surfPts2(:,3),32,pointInd);
% plot3(surfPts2(path,1),surfPts2(path,2),surfPts2(path,3),'r')
% plot3(CL(:,1),CL(:,2),CL(:,3),'b');
% plot3(CLsmoothe(:,1),CLsmoothe(:,2),CLsmoothe(:,3),'g');
% hold off;
% axis equal
% colorbar
% % 

% other outputs that the function expects
% zsq is a binarized image of the outline of the cell
% zbw is a skeletonized version of the cell in xy
%%
zsq = max(hyper_stack,[],3);
zbw = histcn(xyzs,1:size(hyper_stack,1),1:size(hyper_stack,2),[-inf,inf]);
zbw = sum(zbw,3);
zbw = logical(zbw);
zbw = bwmorph(zbw,'skel',inf);
zbw = interp2(double(zbw),upFactor,'*linear'); % BPB is unsure about the mode *linear?

