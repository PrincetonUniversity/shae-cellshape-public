function F=surfaceResample(F,upperLimit,lowerLimit)
% Resamples surface in order to make triangulation more regular. When edge
% lengths are greater than upper limit, the two triangles which are
% connected to that edge are either switched (same quadralateral is made
% but the diagnol which cuts it into two triangles is switched, or the 2
% triangles are cut into 4. 
Fold=F;
Fold2=F;
showflag=0;
for iterations=1:3
    faces2kill=0;
F=faceNormals(F); 
Nvert=size(F.vertices,1);

[E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
%eSelect=find(Elength>l,1,'first');

edges2kill=E(Elength<lowerLimit,:);

while size(edges2kill,1)
%   size(edges2kill,1)
    if showflag
    if exist('h','var') && all(ishandle(h)), delete(h); end
          h = patch(F,'facecolor','w','facealpha',1);  axis equal;  view(3);
       drawnow;
    end
    smallEdgeLength=Elength(Elength<lowerLimit);

V2Fmat=Face2Vert(F.faces); %find connectivity of faces to vert
        iVert=find(smallEdgeLength==min(smallEdgeLength),1,'first'); % index of longest vert

    E1=edges2kill(iVert,1); %identify the two ends
    E2=edges2kill(iVert,2);
    E1=E1(1);E2=E2(1);
    %find the two vertices that are connected to that edge
   if (E1>size(F.vertices,1)) && ...
        (E2>size(F.vertices,1))
overTri=0;
   else
    %find the two vertices that are connected to that edge
   
    overTri=and(V2Fmat(:,E1),V2Fmat(:,E2));
   end
   if sum(overTri(:))==2
    T=F.faces(overTri,:);
Overts=T(T~=E1 & T~=E2);
%find length of other diagnol on the quadrilateral made by the two
%triangles
E1Vert=sqrt(sum((F.vertices(Overts(1),:)-F.vertices(E1,:)).^2,2))+...
    sqrt(sum((F.vertices(Overts(2),:)-F.vertices(E1,:)).^2,2));
E2Vert=sqrt(sum((F.vertices(Overts(1),:)-F.vertices(E2,:)).^2,2))+...
    sqrt(sum((F.vertices(Overts(2),:)-F.vertices(E2,:)).^2,2));


%if both diagnols are longer than the limit, add a vertex in the middle of
%the diagnols, otherwise, switch diagnols

%    F.faces(overTri,:)=[];
    F.faces(F.faces==E2)=E1;
    faces2kill=faces2kill+overTri;
F.vertices(E1,:)=F.vertices(E1,:)/2+F.vertices(E2,:)/2;
       F.faces(faces2kill>0,:)=[];
       faces2kill=0;
   else

       F=faceNormals(F);


    end
[E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
%eSelect=find(Elength>l,1,'first');

edges2kill=E(Elength<lowerLimit & Elength>0,:);
%size(edges2kill,1)
end
F.faces(faces2kill>0,:)=[];
F=faceNormals(F);
[E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
%eSelect=find(Elength>l,1,'first');
valence=sum(vertNeighbours(F.faces,1)>0);
underConnected=any(valence(E)==4,2);
edges2split=E(Elength>upperLimit,:);
%loop through long edges

while size(edges2split,1)
       edges2split= edges2split(1,:);

    if showflag
        if exist('h','var') && all(ishandle(h)), delete(h); end
          h = patch(F,'facecolor','w','facealpha',1);  axis equal;  view(3);
       drawnow;
    end
    edgeLength=Elength(Elength==max(Elength));
    edgeLength=edgeLength(1);
edgeVert1=F.vertices(edges2split(1,1),:);
edgeVert2=F.vertices(edges2split(1,2),:);
newVert=edgeVert1/2+edgeVert2/2;

    V2Fmat=Face2Vert(F.faces); %find connectivity of faces to vert

    E1=edges2split(1); %identify the two ends
    E2=edges2split(2);
    %find the two vertices that are connected to that edge
    overTri=and(V2Fmat(:,E1),V2Fmat(:,E2)); 
    if sum(overTri(:))==2
    T=F.faces(overTri,:);
Overts=T(T~=E1 & T~=E2);
%find length of other diagnol on the quadrilateral made by the two
%triangles
dVert=sqrt(sum((F.vertices(Overts(1),:)-F.vertices(Overts(2),:)).^2,2));

%if both diagnols are longer than the limit, add a vertex in the middle of
%the diagnols, otherwise, switch diagnols
if dVert>upperLimit
    Nvert=size(F.vertices,1);
    F.vertices(Nvert+1,:)=newVert;

T1=T;
T2=T;
T1(T1==E1)=Nvert+1;
T2(T==E2)=Nvert+1;
F.faces(overTri,:)=T1;
F.faces=cat(1,F.faces,T2);
elseif dVert<edgeLength
    T1=cat(1,Overts,E1)';
    T2=cat(1,Overts,E2)';
    F.faces(overTri,:)=[T1;T2];

end

    else
        F=faceNormals(F);
    end

    [E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
%eSelect=find(Elength>l,1,'first');
valence=sum(vertNeighbours(F.faces,1)>0);
underConnected=any(valence(E)==4,2);
edges2split=E(Elength==max(Elength) & Elength>upperLimit,:);

    
end
if size(F.faces,1)==size(Fold.faces,1)
    if all(F.faces==Fold.faces)
        break
    end
end
if size(F.faces,1)==size(Fold2.faces,1)
    if all(F.faces==Fold2.faces)
        break
    end
end

Fold2=Fold;
Fold=F;

end
    
F=faceNormals(F);



