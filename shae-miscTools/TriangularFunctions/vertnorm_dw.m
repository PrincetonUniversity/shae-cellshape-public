function nv=vertnorm_dw(tri,p)
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

allCen=tricentroid(p,tri);

for j=1:length(p)
    % Find centroids and weight based on distance:
    ind01=Vmat(j,:);
    cen=allCen(ind01,:);
    nc=size(cen);
    distcen = cen-repmat(p(j,:),nc(1),1);
    w = 1./sqrt(sum(distcen.*distcen,2));
%     nvecx=mean(w.*nvec(ind01,1));
%     nvecy=mean(w.*nvec(ind01,2));
%     nvecz=mean(w.*nvec(ind01,3));
    % Nomalize
    nvecall=mean(bsxfun(@times,w,nvec(ind01,:)),1);

    lnv = sqrt(sum(nvecall.^2));
    nv(j,:) = nvecall/lnv; 
end
nv=-nv;

end % End vertnorm

function invTRI=buildInverseTriangulation(TRI,neighbourhood)
% Building the inverse triangulation, i.e. a link from node indexes to
% triangle indexes.
%% original
F2Vmat=Face2Vert(TRI);
if nargin==1 ||neighbourhood==1
    neighbourhood=1;
else
    F2Vmat=(F2Vmat*vertNeighbours(TRI,neighbourhood)>0);
    
    
end



[row,col]=find(F2Vmat);

neighbournumber=sum(F2Vmat);

invTRI=zeros(length(neighbournumber),max(neighbournumber));
for i=1:length(neighbournumber)
    invTRI(i,1:neighbournumber(i))=row(col==i);
end


end % End buildInverseTriangulation

%%
function out=removeD0(x)
% Removing duplicate and null values
s=sort(x);
s1=[0,s];
s2=[s,s(length(s))];
ds=(s1-s2);
out=s(logical(ds~=0));
end% End removeD0


%%


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
