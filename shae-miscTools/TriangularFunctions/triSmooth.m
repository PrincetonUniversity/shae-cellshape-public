function Fc2=triSmooth(varargin)
% trismooth takes connectivity between nodes, Tri, and smooths values at
% those nodes. Each point is replaced with the mean of its 1 connected
% neighbours.
% Fc2=triSmooth(F) smooths patch object F, with fields faces and vertices;
%Fc2=trismooth(Tri,Fc) smooths grid with vertices given by tri and faces
%given by Fc, output is JUST VERTICES;
if nargin==1;
    F=varargin{1};
    Tri=F.faces;
    Fc=F.vertices;
elseif nargin==2
    Tri=varargin{1};
    Fc=varargin{2}; 
end


V2Fmat=Face2Vert(Tri)';
V2Vmat=V2Fmat*V2Fmat';
V2Vmat=V2Vmat+speye(size(V2Vmat))>0;

Fc2=Fc;
% BPB comments out to remove smoothing temporarily
for i=1:length(Fc)
    
    Fc2(i,1)=mean(Fc(V2Vmat(i,:)));
end

if nargin==1
    F.vertices=Fc2;
    Fc2=F;
end

    


