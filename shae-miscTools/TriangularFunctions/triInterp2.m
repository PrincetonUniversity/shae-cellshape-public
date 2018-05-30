function [F2,Area2]=triInterp2(F,Area,n)

%interpolates on a close triangular mesh adding points at the center of
%each triangle and redoing connections to maintain a similar connectivity.
%Repeats n times.


vertexNumber=length(F.vertices);
faceNumber=length(F.faces);
E=[F.faces(:,[2,3]);F.faces(:,[3,1]);F.faces(:,[1,2])];
E=sort(E,2);
[E,~,~]=unique(E,'rows');
% Eind=sub2ind([max(E(:)),max(E(:))],E(:,1),E(:,2));
% Emap=zeros(max(E(:)));
% Emap(Eind)=1;
% Emap=(Emap+Emap')>0;
V2Fmat=Face2Vert(F.faces);


fxyz=F.vertices(F.faces,:);
fxyz=reshape(fxyz,[],3,3);

Fsub=squeeze(mean(fxyz,2));
newInd=1:size(Fsub,1);
newInd=newInd'+size(F.vertices,1);

tri1=[newInd,F.faces(:,[2,3])];
tri2=[F.faces(:,1),newInd,F.faces(:,3)];
tri3=[F.faces(:,[1,2]),newInd];
F2=F;
F2.vertices=cat(1,F2.vertices,Fsub);
F2.faces=cat(1,tri1,tri2,tri3);




%after created all new triangles, (3 triangles from a single) split all
%lines from the original mesh by redrawing the connetivity at each old edge
%to go to the 2 closest newly generated points. For now I'm doing a search
%to find the new points, but if I keep track of how I generate each edge,
%this may get rid of the search and be much faster. 
for i=1:length(E)
%T=find(sum(ismember(F2.faces,E(i,:)),2)==2,2,'first');
%T=find(ic==i,2,'first');
%overlap=and(Emap(:,E(i,1)),Emap(:,E(i,2)));
overTri=and(V2Fmat(:,E(i,1)),V2Fmat(:,E(i,2)));
%if sum(overTri)==2
ind=find(overTri);
    
    
    
    
T=F.faces(ind,:);
T=(T~=E(i,1) & T~=(E(i,2)));
ind=mod(find(T')-1,3)*faceNumber+ind;
Ttemp2=F2.faces(ind,:);
newEdge=Ttemp2(Ttemp2~=E(i,1)&Ttemp2~=E(i,2))';
F2.faces(ind,:)=[newEdge,E(i,1);newEdge,E(i,2)];

%end
end


if nargout==2 || nargin>2;
    if nargin==1
[~,Area,~]=triArea(F);
    end
    if nargin>1
        if isempty(Area)
            [~,Area,~]=triArea(F);

        end
        
    end
    
    connectivity=sum(V2Fmat);
    Area2=[Area/3;sum(Area(F.faces)./connectivity(F.faces),2)*2/3];
    
end
%repeat recursivly n times
if nargin ==3
    for i=1:(n-1)
        [F2,Area2]=triInterp2(F2,Area2);
    end
end


