function out=tricurv_v02(tri,p)
% Function to calculate the principal curvatures, and their respective
% directions on a triangular mesh. Approximations of curvature are based on 
% local (N=1) neighborhood elements and vertices.
% Note that calculations at vertices with few adjacent triangles, and hence
% few adjacent vertices, are expanded to a greater neighborhood.
%
%
% Reference:
% 1) Chen and Schmitt (1992) Intrinsic surface properties from surface
% triangulation
% 2) Dong et al. (2005) Curvature estimation on triangular mesh,
%       JZUS
%
% This code makes use of routines: buildInverseTriangualtion.m & removeDO.m
% initially written by: David Gringas. He is gratefully acknowledged
%
% Input:        t <mx3> array of triangle indices into:
%               xyz <nx3> coordinates
% Output:       structure containing: Principal curvatures and directions,
%               Chen and Schmitt's coefficients, etc.
% Version:      1
% JOK 030809

% I/O check:
if nargin~=2
    error('Wrong # of input')
end
if nargout ~= 1
    error('Output is a strcture, wrong designation!')
end
nt=  size(tri);
if nt(2)~=3
    error('Triangle element matrix should be mx3!')
end
mp = size(p);
if mp(2)~=3
    error('Vertices should be nx3!')
end


% 1) Compute vertex normals (weighted by distance)
nv=vertnorm_dw(tri,p);

% 2) Define neighborhoods and calculate weighted normals
% Average at vertices
itri=buildInverseTriangulation(tri,1);

for j=1:length(p)
    ind01=removeD0(itri(j,:));% Triangles adjacent to vertex j
    ind02 = tri(ind01,:);ind02 = unique(ind02(:));% Vertices of adjacent triangles
    if length(ind02)<6 % Larger neighborhood for triangles with less than 5 unique triangles
        for i=1:length(ind02)
            indaux01 = removeD0(itri(ind02(i),:));% Triangles adjacent to neighborhood vertices in ind02(i)
            indaux02 = tri(indaux01,:);% Vertices  of larger neighborhood
            indaux02 = unique(indaux02(:));
            ind02 = unique([ind02;indaux02]);
        end
    end
    % Define tangent vectors
    pdist = p(ind02,:)-repmat(p(j,:),length(ind02),1); %displacement vectors to neighbours
    aux01 = pdist-repmat(dot(pdist,nv(ind02,:),2),1,3).*nv(ind02,:);%project away normal component
    
    t = aux01./repmat(sqrt(sum(aux01.*aux01,2)),1,3);%normalize
    % Normal curvature
    kn = -dot(pdist,nv(ind02,:)- ...
                repmat(nv(j,:),length(ind02),1),2)./dot(pdist,pdist,2);
    [mkn,ind03] = max(kn);
    % Local principal directions CHANGE THIS TO FIT ELLIPSE TO KN
    out.e1(j,:) = t(ind03,:);
    aux02 = cross(out.e1(j,:),nv(j,:),2);
    out.e2(j,:) = aux02./repmat(sqrt(sum(aux02.*aux02,2)),1,3);
    % Chen and Schmitt: calculation of coefficients
    thetai = real(acos(dot(t,repmat(out.e1(j,:),length(ind02),1),2)));
    % Check for NaN
    thetai(isnan(thetai)) = 0; % Can be set to zero as all coefficients 
    % contain a sin term and thus are zero
    
    c2=cos(thetai).^2;
    s2=sin(thetai).^2;
    cs=cos(thetai).*sin(thetai);
    
    
    
    cen=p(ind02,:);
    nc=size(cen);
    distcen = cen-repmat(p(j,:),nc(1),1);
    w = sqrt(sum(distcen.*distcen,2));
    w(w==0)=Inf;
%     w(w<3)=3;
%     w=exp(-(w/2.5).^2/2);
   w=1./w;
    w(isinf(w))=0;
    W=diag(w);
    %W=W>0;
    A=[c2,cs,s2];
    Ainv=(pinv(A)*W*A)\pinv(A)*W; %weighted pseudo inverse matrix
   % Ainv=pinv(A);
    
    
    kn(isnan(kn))=0;
    
    abc=Ainv*kn;
    a(j,1)=abc(1);
    b(j,1)=abc(2);
    c(j,1)=abc(3);
    
%     if any(any([.5*(a+c+sqrt((a-c).^2+4*b.^2)),.5*(a+c-sqrt((a-c).^2+4*b.^2))]...
%             <-.4))
%         keyboard
%     end


end

% Gaussian, mean normal, principal curvatures
aux03 = sort([.5*(a+c+sqrt((a-c).^2+4*b.^2)),.5*(a+c-sqrt((a-c).^2+4*b.^2))],2);
out.k1 = aux03(:,1);
out.k2 = aux03(:,2);
out.km = .5*(out.k1+out.k2);
if mean(out.km)<0
    out.km=-out.km;
    out.k1=-aux03(:,2);
    out.k2=-aux03(:,1);
end

out.kg = out.k1.*out.k2;

end % End tricurv_v01




%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
function invTRI=buildInverseTriangulation(TRI,neighbourhood)
% Building the inverse triangulation, i.e. a link from node indexes to
% triangle indexes.
%% original
F2Vmat=Face2Vert(TRI);
if nargin==1 ||neighbourhood==1
    neighbourhood=1;
else
    F2Vmat=(F2Vmat*vertNeighbours(TRI,neighbourhood))>0;
    
    
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
itri=buildInverseTriangulation(tri,3);
nvecx=zeros(length(p),1);
nvecy=zeros(length(p),1);
nvecz=zeros(length(p),1);
for j=1:length(p)
    % Find centroids and weight based on distance:
    ind01=removeD0(itri(j,:));
    cen=tricentroid(p,tri(ind01,:));
    nc=size(cen);
    distcen = cen-repmat(p(j,:),nc(1),1);
    w = 1./sqrt(sum(distcen.*distcen,2));
    nvecx=mean(w.*nvec(ind01,1));
    nvecy=mean(w.*nvec(ind01,2));
    nvecz=mean(w.*nvec(ind01,3));
    % Nomalize
    lnv = sqrt(nvecx^2+nvecy^2+nvecz^2);
    nv(j,:) = [nvecx/lnv,nvecy/lnv,nvecz/lnv]; 
end

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


%%
% function Vmat=vertNeighbours(TRI,neighbourhood)
% 
% if nargin==1
% neighbourhood=1;
% end
% vertexNumber=max(TRI(:));
% Vmat=zeros(vertexNumber);
% for i=1:vertexNumber
% [row,col]=find(TRI==i);
% neighpoints=TRI(row,:);
% neighpoints(neighpoints==i)=[];
% neighpoints=unique(neighpoints);
% Vmat(neighpoints,i)=1;
% end
% 
% Vmat=Vmat^neighbourhood;
% Vmat=Vmat.*(1-eye(size(Vmat)));
% Vmat=Vmat>0;
% % [row,col]=find(Vmat);
% % 
% % neighbournumber=sum(Vmat);
% % 
% % invTri=zeros(vertexNumber,max(neighbournumber));
% % for i=1:vertexNumber;
% %     invTri(i,1:neighbournumber(i))=row(col==i);
% % end
% % end
% end

%%


function Fmat=faceNeighbours(TRI,neighbourhood)
if nargin==1
neighbourhood=1;
end

vertexNumber=max(TRI(:));
faceNumber=length(TRI);


Fmat=zeros(faceNumber);
for i=1:faceNumber
    neighpoints=TRI(i,:);
[neighbourTri,col]=find(ismember(TRI,neighpoints));
neighbourTri(neighbourTri==i)=[];
neighbourTri=unique(neighbourTri);
Fmat(neighbourTri,i)=1;
end


Fmat=Fmat^neighbourhood;
end

% %%
% function F2Vmat=Face2Vert(TRI)
% vertexNumber=max(TRI(:));
% faceNumber=length(TRI);
% 
% F2Vmat=zeros(faceNumber,vertexNumber);
% for i=1:faceNumber
%     neighpoints=TRI(i,:);
% F2Vmat(i,neighpoints)=1;
% end
% 
% 
% end