function Vmat=vertDistMat(F,path,option)
%This function takes a mesh given by the patch object F finds the distance
%to its nearest neighbours. The output of Vmat is a matrix in which
%Vmat(i,j) is the distance between vertex i and j if they are connected. If
%they are not connected  Vmat(i,j)=0. 

%there is also an option to specify a path and an option. This will
%constrain distances in a certain way, either to go in the direction of a curve 
% or circle around it. Path is a matrix [x,y,z] of
%coordinates specifying a curve. Distances in Vmat can be constrained to
%either by in the direction of the path ('linear') or perpendicular to it
%('circular'). 

%if path is specified, it 
if nargin==1
    option='none';
end


vertexNumber=max(F.faces(:));

E=[F.faces(:,[2,3]);F.faces(:,[3,1]);F.faces(:,[1,2])];
[E,~,~]=unique(E,'rows');

E=sort(E,2);
E=[E;E(:,2),E(:,1)];

P1=F.vertices(E(:,1),:);
P2=F.vertices(E(:,2),:);
Dvec=P1-P2;
D=sqrt(sum(Dvec.^2,2));

if ~strcmp(option,'none')
    
[~,path_vectors]=gradient(path);
[~,ind]=pdist2(path,F.vertices,'Euclidean','Smallest',1);
Ovec=F.vertices-path(ind,:);
Ovec=Ovec(E(:,1),:);
Tvec=path_vectors(ind,:);
Tvec=Tvec(E(:,1),:);
end

switch option
    case 'circular'
        s=cross(Ovec,Tvec);
        s=normalizeRows(s);
        Dvec=normalizeRows(Dvec);
        theta=acos(dot(Dvec, s,2));
        S=theta<80*pi/180;
        D=D.*S;
        D=double(D);

        Vmat=sparse(E(:,1),E(:,2),D,...
    vertexNumber,vertexNumber);
    case 'linear'
        D=D.*(sign(dot(Dvec,Tvec))>-1);
        D(D<0)=0;
        D=double(D);

        Vmat=sparse(E(:,1),E(:,2),D,...
    vertexNumber,vertexNumber);
       otherwise
        D=double(D);

        Vmat=sparse(E(:,1),E(:,2),D,...
    vertexNumber,vertexNumber);
end
        


% 
% Vmat=sparse([E(:,1),E(:,2)],[E(:,2),E(:,1)],[D,D],...
%     vertexNumber,vertexNumber);

% 
% for i=1:length(TRI);%vertexNumber
%     
%     
% % [row,col]=find(TRI==i);
% % neighpoints=TRI(row,:);
% % neighpoints(neighpoints==i)=[];
% % neighpoints=unique(neighpoints);
% 
% 
% V1=TRI(i,:);
% V2=V1([3,1,2]);
% neighCoord=VERT(V1,:);
% neighDist=neighCoord-neighCoord([3,2,1],:);
% neighDist=sqrt(sum(neighDist.^2,2));
% Vrow=[V1,V2];
% Vcol=[V2,V1];
% Vind=vertexNumber*(Vcol-1)+Vrow;
% Vmat(Vind)=[neighDist,neighDist];
% 
% 
% 
% 
% 
% 
% 
% end
% 
