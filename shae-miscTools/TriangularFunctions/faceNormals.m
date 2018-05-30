function [F,Vout]=faceNormals(F)


%function takes a patch object and computs the normal to each face using
%the dot products of all the vectors.

%%remove duplicate points in same triangle
F.faces(any(F.faces==circshift(F.faces,[0,1]) | F.faces==circshift(F.faces,[0,2]),2),:)...
    =[];

%% remove duplicate triangles

F.faces=sort(F.faces,2);
F.faces=unique(F.faces,'rows');
%% remove edges that only neighbour 1 triangle until no triangle has this problem
H=1;

while  H~=0
TRI=F.faces;
E=[TRI(:,[1,3]);TRI(:,2:3);TRI(:,1:2)];
E=sort(E,2);
[~,ia,ib]=unique(E,'rows');
H=histc(ib,1:1:max(ib(:)));

underConnected=ia(H<2);
[triKill,~]=ind2sub([length(TRI),3],underConnected);
overConnected=ia(H>2);

for i=1:length(overConnected);
    F.vertices(E(overConnected(i),1),:)=F.vertices(E(overConnected(i),1),:)/2+...
    F.vertices(E(overConnected(i),2),:)/2;
    F.faces(F.faces==E(overConnected(i),2))=E(overConnected(i),1);
    triKill=cat(1,triKill,find(sum(F.faces==E(overConnected(i),1),2)==2));

end
F.faces(triKill,:)=[];
H=numel(overConnected)+numel(underConnected);
end

%% remove vertices that are not connected to other vertices
Vmat=vertNeighbours(F.faces,1);
Vout=(sum(Vmat,2)<=1);
change=cumsum(Vout);
%change(indOut)=[];
F.vertices(max(F.faces(:))+1:size(F.vertices,1),:)=[];
F.vertices(Vout,:)=[];
F.faces=F.faces-change(F.faces);

%% Orient faces in all the same direction
faceNum=length(F.faces);
FedgeMap=faceNeighbours(F.faces,1,1);
T0=zeros(faceNum,1);
oldTri=T0;   %faces we've looked at. start with just one.
T0(1)=1;
Tri1=F.faces(1,:);
tic
while ~all(oldTri) && toc<10; %go through all faces

T1=FedgeMap*T0;   %look at neighbours of points from previous round
T1=T1.*~oldTri;     %excluding points we've already looked at
Tri2=F.faces(T1>0,:);
%triangles that need switching will be equal to one of the 3 rotated
%permutations of the original indices
    Tri11=cat(1,Tri1,Tri1(:,[2,3,1]),Tri1(:,[3,1,2])); 
    Tri11=permute(Tri11,[3,2,1]); %turn into 3d matrix for bsxfun
% if they have 2 items equal, that means they share an edges 
% in the same orientation and should be switched
T=(bsxfun(@eq,Tri11,Tri2));
T=any(sum(T,2)==2,3);

    T=any(T,3); % triangles need switching
    Tri2Flip=Tri2(:,[2,1,3]);
    Tri2(T,:)=Tri2Flip(T,:);
%reinsert changes into original triangulation
    F.faces(T1>0,:)=Tri2;
Tri1=Tri2;
oldTri=(T1+oldTri)>0;
T0=T1;
end

%make sure faces are oriented away from center of mass when crossed in the
%order 12 cross 23 or 23 cross 31 or 12 cross 13
p=F.vertices;
tri=F.faces;
normvec=vertnorm(tri,p);
p=bsxfun(@minus,p,mean(p,1));
p=normalizeRows(p);
normvec=normalizeRows(normvec);



 if sum(dot(p,normvec,2))<0
    F.faces=F.faces(:,[2,1,3]);
%     display('switch');
 end
    
Vmat=vertNeighbours(F.faces,1);
indOut=(sum(Vmat,2)<=1);
change=cumsum(indOut);
%change(indOut)=[];
F.vertices(max(F.faces(:))+1:size(F.vertices,1),:)=[];
F.vertices(indOut,:)=[];
F.faces=F.faces-change(F.faces);

