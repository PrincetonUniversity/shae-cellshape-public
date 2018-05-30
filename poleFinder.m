function [poleIdx,isPolarRegionVertex,isPolarRegionFace] = poleFinder(F,Fcoord,prefs);
% this function will attempt to find, geometrically, regions that would be
% considered the pole
% [poleIdx,isPolarRegionVertex,isPolarRegionFace] = poleFinder(F,Fcoord,prefs);
if (nargin<2)
    
    F = [];
end

if nargin <3
    prefs = [];
end

if isempty(F)
    % without any inputs, prompt the user to select a file
    [filename, filePath] = uigetfile();
    load(fullfile(filePath,filename))
end

if (not(isstruct(F)))
    % establish some defaults
    
    % 1 displays a new figure with poles, 0 doesn't display anything
    prefs.displayFigure = 1;
    
end


% establish default values for fields that are missing
if not(isfield(prefs,'radiusMultiplier'))
    prefs.radiusMultiplier = 1.5;
end

if not(isfield(prefs,'displayFigure'))
    prefs.displayFigure = 1;
end

if not(isfield(prefs,'flipChirality'))
    prefs.flipChirality = -1;
end

if not(isfield(prefs,'cdataName'))
    prefs.cdataName = 'gc';
end


%%
% find the index that is at the pole
poleCoord1 = Fcoord.CL(1,:);
distanceToPole = bsxfun(@minus,F.vertices,poleCoord1);
distanceToPole = sqrt(sum(distanceToPole.^2,2));
[~,pnt1] = min(distanceToPole);
%

% look near each of these and build a neighborhood/patch that is below a
% certain radius
[F2Vmat,invTri]=Face2Vert(F.faces);
[~,Amesh] = triArea(F);
[r0,sigma,excessProb] = fitCellRadius(Fcoord.r(:),Amesh(:));

% distance along the centerline that these surface elements map to
Ldist = Fcoord.L(:);

% maximum distance from pole along surface to include in the pole
maxRadiusInPole = r0*prefs.radiusMultiplier;

% start at pnt1
isAnotherRound = true;
% create 'sparse' array to define surface region
currentPatch = sparse(size(F2Vmat,2),1);
% start at a point
currentPatch(pnt1) = true;
% vertex to vertex connectivity
Vert2Vert = F2Vmat'*F2Vmat;
Vert2Vert(Vert2Vert>0) = 1;
% only keep vertices that are close to the pole
threshRadius = Ldist<(maxRadiusInPole);
patchSize = 1; % counts vertices
while isAnotherRound
    % increase size of patch to include one more set of neighbors
    currentPatch = logical((currentPatch'*Vert2Vert).*threshRadius')';
    % check to see if there were any new neighbors added
    isAnotherRound = nnz(currentPatch)-patchSize;
    patchSize = nnz(currentPatch);
end

% look up which faces are included in this new surface patch
facesKept = invTri(currentPatch,:);
facesKept = unique(facesKept);
facesKept(facesKept==0) = [];


%% start at pnt2
% find the index that is at the pole
poleCoord2 = Fcoord.CL(end,:);
distanceToPole = bsxfun(@minus,F.vertices,poleCoord2);
distanceToPole = sqrt(sum(distanceToPole.^2,2));
[~,pnt2] = min(distanceToPole);
%

% look near each of these and build a neighborhood/patch that is below a
% certain radius

% distance along the centerline that these surface elements map to
Ldist = abs( Fcoord.L(pnt2)-Fcoord.L(:));


% start at pnt2
isAnotherRound = true;
% create 'sparse' array to define surface region
currentPatch = sparse(size(F2Vmat,2),1);
% start at a point
currentPatch(pnt2) = true;
% vertex to vertex connectivity
Vert2Vert = F2Vmat'*F2Vmat;
Vert2Vert(Vert2Vert>0) = 1;
% only keep vertices that are close to the pole
threshRadius = Ldist<(maxRadiusInPole);
patchSize = 1; % counts vertices
while isAnotherRound
    % increase size of patch to include one more set of neighbors
    currentPatch = logical((currentPatch'*Vert2Vert).*threshRadius')';
    % check to see if there were any new neighbors added
    isAnotherRound = nnz(currentPatch)-patchSize;
    patchSize = nnz(currentPatch);
end

% look up which faces are included in this new surface patch
facesKept2 = invTri(currentPatch,:);
facesKept2 = unique(facesKept2);
facesKept2(facesKept2==0) = [];

% build a new patch style object with these
F3 = F;
    F3.vertices(:,2) = prefs.flipChirality*F3.vertices(:,2);

F4 = F;
    F4.vertices(:,2) = prefs.flipChirality*F4.vertices(:,2);

F4.faces = F4.faces(cat(1,facesKept2,facesKept),:);
isPolarRegionFace = ismember(F3.faces,F4.faces,'rows');
F3.faces(isPolarRegionFace,:) = [];


%% create image
if prefs.displayFigure
    gcf;
    clf;
    set(gcf,'Color','w');
    Ftemp = F;
    Ftemp.vertices(:,2) = prefs.flipChirality*Ftemp.vertices(:,2);
    
    % choose what to color the surface based on
    switch prefs.cdataName
        
        case {'gc','gaussian'}
            Fcurvature = tricurv_v02(Ftemp.faces,Ftemp.vertices);
            Fcurvature.kg=triSmooth(F.faces,Fcurvature.kg);

            cdata = Fcurvature.kg(:);
        case {'r','radius'}
            cdata = Fcoord.r(:);
    end
    
    
    % display this resulting patch and pole
    patch(F3,'facevertexCdata',cdata,'facecolor','interp','facealpha',0.2,'edgecolor','none');
    
    axis equal;
    axis vis3d;
    hold on;
    
    patch(F4,'facevertexCdata',cdata,'facecolor','interp','facealpha',0.8,'edgecolor','none');
    scatter3(Ftemp.vertices(pnt2,1),Ftemp.vertices(pnt2,2),Ftemp.vertices(pnt2,3),250,'k*');
    scatter3(Ftemp.vertices(pnt1,1),Ftemp.vertices(pnt1,2),Ftemp.vertices(pnt1,3),250,'k*');
end

% in invTri, the row index is the name of a vertex and values in the
% columns are the faces that that vertex touches

% vertices in the polar region are where ANY of the faces in the polar
% region touch them
invTri(invTri==0) = nan;
[invTriRows,invTriCols] = size(invTri);

invTri2 = invTri(:);
isPolarRegionVertex = ismember(invTri2,find(isPolarRegionFace));

isPolarRegionVertex = reshape(isPolarRegionVertex,invTriRows,invTriCols);
isPolarRegionVertex = any(isPolarRegionVertex,2);

% scatter3(Ftemp.vertices(isPolarRegionVertex,1),Ftemp.vertices(isPolarRegionVertex,2),Ftemp.vertices(isPolarRegionVertex,3),250,'ro');

poleIdx = [pnt1,pnt2];




