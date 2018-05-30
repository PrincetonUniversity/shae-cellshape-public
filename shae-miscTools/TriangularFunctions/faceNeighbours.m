function Fmat=faceNeighbours(TRI,neighbourhood,edgeflag)
switch nargin
    case 1
neighbourhood=1;
edgeflag=0;
    case 2
        edgeflag=0;
end

if edgeflag
    faceNum=length(TRI);
E=[TRI(:,[2,3]);TRI(:,[3,1]);TRI(:,[1,2])];
E=sort(E,2);

[~,ia,ic]=unique(E,'rows','first');
[~,ia2,ic2]=unique(E,'rows','last');
ia=mod(ia,faceNum);
ia(ia==0)=faceNum;
ia2=mod(ia2,faceNum);
ia2(ia2==0)=faceNum;
Fmat=sparse(ia,ia2,ones(size(ia2)),faceNum,faceNum);
Fmat=Fmat+Fmat';

else


Fmat=Face2Vert(TRI);
Fmat=Fmat*Fmat';


end


Fmat=Fmat^neighbourhood>0;
