function nv=vertnorm(tri,p)
% Function to compute normal vector of vertices comprising a triangular
% mesh. Based on trinormal and computeNormalVectorTriangulation.m by David
% Gringas
% Input:    tri mx3 <triangular index matrix>
%           p nx3 array of vertices
% Output:   nvec <nx3> array of normal vectors
% JOK230709
% Version: 1

% I/O check
% to be completed

% Construct vectors 
v = [p(tri(:,3),1)-p(tri(:,1),1), p(tri(:,3),2)-p(tri(:,1),2), p(tri(:,3),3)-p(tri(:,1),3)];
w = [p(tri(:,2),1)-p(tri(:,1),1), p(tri(:,2),2)-p(tri(:,1),2), p(tri(:,2),3)-p(tri(:,1),3)];


% Calculate cross product
normvec = [v(:,2).*w(:,3)-v(:,3).*w(:,2), ...
    -(v(:,1).*w(:,3)-v(:,3).*w(:,1)), ...
    v(:,1).*w(:,2)-v(:,2).*w(:,1)];
% Normalize
lnvec = sqrt(sum(normvec.*normvec,2));
nvec = normvec./repmat(lnvec,1,3);

% Average at vertices
%itri=buildInverseTriangulation(tri,1);
Vmat=Face2Vert(tri)'>0;
D=speye(length(p));
LL=Vmat'*D;
allCen=tricentroid(p,tri);

allX=bsxfun(@times,LL,normvec(:,1));
allY=bsxfun(@times,LL,normvec(:,2));
allZ=bsxfun(@times,LL,normvec(:,3));


nv=[sum(allX)',sum(allY)',sum(allZ)'];
nv=-normalizeRows(full(nv));

end % End vertnorm



%%
function out = tricentroid(v,tri)
% Function to output the centroid of triangluar elements.
% Note that the output will be of length(tri)x3
% Input:    <v>     nx2 or 3: vertices referenced in tri
%           <tri>   mx3: triangle indices
% Version:      1
% JOK 300509

% I/O check
[nv,mv]=size(v);
[nt,mt]=size(tri);
if mv==2
    v(:,3) = zeros(nv,1);
elseif mt~=3
   tri=tri';
end

out(:,1) = 1/3*(v(tri(:,1),1)+v(tri(:,2),1)+v(tri(:,3),1));
out(:,2) = 1/3*(v(tri(:,1),2)+v(tri(:,2),2)+v(tri(:,3),2));
out(:,3) = 1/3*(v(tri(:,1),3)+v(tri(:,2),3)+v(tri(:,3),3));
end%tricentroid
