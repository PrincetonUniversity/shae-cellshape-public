function [cellCoordOut, rotMat] = centerShape(shapeCoord)
% CENTERSHAPE rotates and translates to [0,0,0] and aligns to [x,y,z]
%
% [c] = centerShape(shapeCoord) returns the recentered object.
% 
% [c,r] = centerShape(shapeCoord) returns the recentered object and the
% rotation used.
%
% Example:
%   [xx,yy,zz] = sphere(22);
%   shapeCoord = [xx(:),yy(:),zz(:)];
%   shapeCoord = unique(shapeCoord,'rows');
%   shapeCoord = [shapeCoord(:,1),shapeCoord(:,2)*0.5,shapeCoord(:,3)*5];
%   [c] = centerShape(shapeCoord);
%   clf;
%   scatter3(c(:,1),c(:,2),c(:,3),'rx');
%   hold on;
%   scatter3(shapeCoord(:,1),shapeCoord(:,2),shapeCoord(:,3),'bo');
%   axis equal;
%
% See also EIGS 

%% find triangulated surface for the test

% [xx,yy,zz] = sphere(22);
% shapeCoord = [xx(:),yy(:),zz(:)];
% shapeCoord = unique(shapeCoord,'rows');
% shapeCoord = [shapeCoord(:,1),shapeCoord(:,2)*2,shapeCoord(:,3)*5];
% shapeCoord = shapeCoord+0.2*randn(size(shapeCoord));
try
dt1 = delaunayTriangulation(shapeCoord);
catch
    dt1 = DelaunayTri(shapeCoord);
end

dt2 = freeBoundary(dt1);

storedAreas1 = zeros(size(shapeCoord,1),1);
x1 = shapeCoord(dt2(:,1),:);
x2 = shapeCoord(dt2(:,2),:);
x3 = shapeCoord(dt2(:,3),:);

areas1 = sqrt(sum(((cross(x3-x2,x3-x1)/2)).^2,2));
areas1 = repmat(areas1,3,1); % each triangle gets added to three vertices
for iVertex = 1:size(shapeCoord,1)
    % for iVertex = 1:10
    storedAreas1(iVertex) = sum(areas1(dt2(:)==iVertex));
end

% fix center of mass

com2x = bsxfun(@times,shapeCoord,storedAreas1);
com2 = sum(com2x,1)./sum(storedAreas1);
shapeCoord = bsxfun(@minus,shapeCoord,com2);


%% calculate mass moments / inertial tensor
for ii=1:3
    for jj=1:3
        inertTen(ii,jj) = sum(sum(shapeCoord(:,ii).*shapeCoord(:,jj).*storedAreas1))/sum(storedAreas1);
    end
end

%% apply rotation
[rotMat, massEig] = eigs(inertTen);
cellCoordOut = shapeCoord*(rotMat);