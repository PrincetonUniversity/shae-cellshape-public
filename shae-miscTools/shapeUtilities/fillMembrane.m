function membrane=fillMembrane(membrane,varargin)


%%%% fills membrane mask given by set of coordinates
% membrane=fillMembrane(membrane,F) 
% or membrane=fillMembrane(membrane,F,upFactor) 

%       works for 3D membrane image and patch object F
% membrane=fillMembrane(membrane,x,y,z)
%or membrane=fillMembrane(membrane,x,y,z,upFactor)
%      works for 3Dmembrane image and x y z matricies for points.
[row,col,stack]=size(membrane);


if nargin<=3
    if nargin==2
        upFactor=1;
    else
        upFactor=varargin{2};
    end
    
   F=varargin{1};
    membrane=membrane>0;
        M2=membrane>0;

    membrane=padarray(membrane,[0,0,1],true,'both');
   
    for dfactor=3:3:min(size(membrane))/2
        
    membrane2=imdilate(membrane,true(dfactor,dfactor,dfactor));
    membrane2=imfill(double(membrane2~=0),'holes');
    membrane2=imerode(membrane2,true(dfactor,dfactor,dfactor));

mdifference=sum((membrane2(:)>0)-membrane(:)>0);
if mdifference>sum(M2(:))/2
    membrane=membrane2;
    break
end
    end
    membrane=membrane(:,:,2:end-1);
    
    
    [mx,my,mz] = ind2sub(size(M2),find(M2));
    membranePix=[mx,my,mz]/upFactor;
    [~,I] = pdist2(F.vertices,membranePix,'euclidean','Smallest',1);
            Nvec=vertnorm_dw(F.faces,F.vertices);
dM=membranePix-F.vertices(I,:);
dVec=Nvec(I,:);
memVal=2*(dot(dM,dVec,2)+1/4);
memVal(memVal<0)=0;
memVal(memVal>1)=1;
memVal(isnan(memVal))=0;
membrane(M2)=memVal;
elseif nargin<=5
    
    if nargin==5
        upFactor=varargin{4};
    else
        upFactor=1;
    end
    
    
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    
    [~,~,V]=surfArea(x,y,z);
    
    V=V;
    F=[reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
    
    
    
    
    V=reshape(V,[],3);
    M2=membrane>0;
    for dfactor=1:3:min(size(membrane))/2
        
    membrane2=imdilate(membrane,true(dfactor,dfactor,dfactor));
    membrane2=imfill(double(membrane2~=0),'holes');
    membrane2=imerode(membrane2,true(dfactor,dfactor,dfactor));

mdifference=sum((membrane2(:)>0)-(membrane(:)>0));
if mdifference>sum(membrane(:)>0)/2
    membrane=membrane2;
    break
end
    end
    
    objects=(membrane+M2==1);
    objectslabel=bwlabeln(objects,6);
    maxVal=0;
    for iObject=1:max(unique(objectslabel));
        if sum(objectslabel(:)==iObject)>maxVal
            maxVal=sum(objectslabel(:)==iObject);
            maxIdx=iObject;
        end
        
    end
    membrane(objectslabel~=maxIdx)=0;
    
    
    [mx,my,mz] = ind2sub(size(M2),find(M2));
    membranePix=[mx,my,mz]/upFactor;
    [~,I] = pdist2(F,membranePix,'euclidean','Smallest',1);
    dM=membranePix-F(I,:);
    V=V(I,:);
    memVal=2*(dot(dM,-V,2)+1/4);
    memVal(memVal<0)=0;
    memVal(memVal>1)=1;
    memVal(isnan(memVal))=0;
    membrane(M2)=memVal;
    
    
end
