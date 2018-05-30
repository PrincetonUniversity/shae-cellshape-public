function [F,VertArea]=surfaceResample4(F,upperLimit,lowerLimit,VertArea)
% Resamples surface in order to make triangulation more regular. When edge
% lengths are greater than upper limit, the two triangles which are
% connected to that edge are either switched (same quadralateral is made
% but the diagnol which cuts it into two triangles is switched, or the 2
% triangles are cut into 4. same as 1 but adds area tracking
Fold=F;
Fold2=F;
showflag=0;
for iterations=1:150
    faces2kill=0;
    F=faceNormals(F);
    Nvert=size(F.vertices,1);
    
    [E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
    %eSelect=find(Elength>l,1,'first');
    
    edges2kill=E(Elength<lowerLimit,:);
    
    
    while size(edges2kill,1)
        if showflag
            if exist('h','var') && all(ishandle(h)), delete(h); end
            h = patch(F,'facecolor','w','facealpha',1);  axis equal;  view(3);
            drawnow;
        end
        smallEdgeLength=Elength(Elength<lowerLimit & Elength>0);
        
        V2Fmat=Face2Vert(F.faces); %find connectivity of faces to vert
        iVert=smallEdgeLength==min(smallEdgeLength); % index of longest vert
        F2Fmat=double(V2Fmat*V2Fmat'>0);
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
            F.vertices(E1,:)=F.vertices(E1,:)/2+F.vertices(E2,:)/2;
            
            F.faces(F.faces==E2)=E1;
            faces2kill=faces2kill+overTri;
            
            VertArea(E1)=VertArea(E1)+VertArea(E2);
            VertArea(E2)=0;
            
            
        else
                F.faces(faces2kill>0,:)=[];
            F=faceNormals(F);
            faces2kill=zeros(length(F.faces),1);
        end
        [E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
        %eSelect=find(Elength>l,1,'first');
        edges2kill=E(Elength<lowerLimit & Elength>0,:);
        %size(edges2kill,1)
    end
    F.faces(faces2kill>0,:)=[];
    
    [F,Vout]=faceNormals(F);
    if size(Vout,1)~=size(VertArea,1)
        VertArea(Vout)=0;
    end
    
    VertArea(VertArea==0)=[];
    
    
    [E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
    %eSelect=find(Elength>l,1,'first');
    valence=sum(vertNeighbours(F.faces,1)>0);
    underConnected=any(valence(E)==4,2);
    
    edges2split=E(Elength>upperLimit,:);
    %loop through long edges
    edgeLength=Elength(Elength>upperLimit);
    
    Nvert=size(F.vertices,1);
    
    while size(edges2split,1)
        
        if showflag
            if exist('h','var') && all(ishandle(h)), delete(h); end
            h = patch(F,'facecolor','w','facealpha',1);  axis equal;  view(3);
            drawnow;
        end
        
        edgeVert1=F.vertices(edges2split(:,1),:);
        edgeVert2=F.vertices(edges2split(:,2),:);
        newVert=edgeVert1/2+edgeVert2/2;
        
        V2Fmat=Face2Vert(F.faces); %find connectivity of faces to vert
        iVert=find(edgeLength==max(edgeLength),1,'first'); % index of longest vert
        
        E1=edges2split(iVert,1); %identify the two ends
        E2=edges2split(iVert,2);
        %find the two vertices that are connected to that edge
        overTri=and(V2Fmat(:,E1),V2Fmat(:,E2));
        if sum(overTri)==2
            T=F.faces(overTri,:);
            Overts=T(T~=E1 & T~=E2);
            %find length of other diagonal on the quadrilateral made by the two
            %triangles
            dVert=sqrt(sum((F.vertices(Overts(1),:)-F.vertices(Overts(2),:)).^2,2));
            
            %if both diagonal are longer than the limit, add a vertex in the middle of
            %the diagonal, otherwise, switch diagonal
            if dVert>upperLimit
                
                F.vertices(Nvert+1,:)=newVert(iVert,:);
                VertArea(Nvert+1)=1/4*sum(VertArea([E1;E2;Overts]));
                VertArea([E1;E2;Overts])=3/4*VertArea([E1;E2;Overts]);
                Nvert=Nvert+1;
                
                T1=T;
                T2=T;
                T1(T1==E1)=Nvert;
                T2(T==E2)=Nvert;
                F.faces(overTri,:)=T1;
                F.faces=cat(1,F.faces,T2);
                
                
            elseif dVert<edgeLength(iVert)
                T1=cat(1,Overts,E1)';
                T2=cat(1,Overts,E2)';
                F.faces(overTri,:)=[T1;T2];
            end
            
            Vmat=vertNeighbours(F.faces,1);
            V2Fmat=Face2Vert(F.faces);
            valence=sum(Vmat>0);
            underConnected=valence==4;
            Vmat2=(Vmat^2>1)-speye(size(Vmat))-Vmat;
            nearUnderC=Vmat2(:,end).*underConnected';
            nearUnderC=(nearUnderC>0).*valence';
            
            if any(nearUnderC);
                E1=find(nearUnderC,1,'first');
                E2=length(F.vertices);
                overTri=and(Vmat(:,E1),Vmat(:,E2));
                overTri=V2Fmat*overTri==2;
                T=F.faces(overTri,:);
                T=T(1,:);
                Overts=T(T~=E1 & T~=E2)';
                T1=cat(1,Overts(1),E1,E2)';
                T2=cat(1,Overts(2),E2,E1)';
                F.faces(overTri,:)=[T1;T2];
            end
            
            
            [E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
            %eSelect=find(Elength>l,1,'first');
            select=(Elength>upperLimit);
            if any(nearUnderC)
                select=and(select,~any((E==E1 |E==E2),2));
            end
            edges2split=E(select,:);
            edgeLength=Elength(select);
            
            
            
        else
            F=faceNormals(F);
            [E,Elength]=EdgeLength(F); %detect pairs of vertices for edges and lengths
            select=(Elength>upperLimit);
            edges2split=E(select,:);
            edgeLength=Elength(select);
            
            
        end
        
        
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




