function Image=coord2image3dfill(varargin)
%Create 3d image mask from coordinates, Coordinates can be in 3 vectors, x
%y and z, 3 matrices, x, y and z, or as a patch object, with vertices and
%faces.
% 
% Image=coord2image3d(F,imsize,upfactor,nm_per_pixel,Fweight);
% Image=coord2image3d(x,y,z,imsize,upfactor,nm_per_pixel);
% upfactor is the factor which the image is scaled up by, nm_per_pixel is
% the amount the coordinates are scaled down by. for x y z as a matrix,
% points are weighted by area. 

if isfield(varargin{1},'vertices')
    F=varargin{1};
    x=F.vertices(:,1);
    y=F.vertices(:,2);
    z=F.vertices(:,3);
    
    upFactor=1;
    nm_per_pixel=1;
    aScale=ones(size(F.vertices,1),1);
    switch nargin
        case 2
            imsize=varargin{2};
            
            
        case 3
            upFactor=varargin{3};
            
            imsize=varargin{2};
        case 4
            imsize=varargin{2};
            upFactor=varargin{3};
            nm_per_pixel=varargin{4};
            
        case 5
            imsize=varargin{2};
            aScale=varargin{5};
            upFactor=varargin{3};
            nm_per_pixel=varargin{4};

    end
else
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    upFactor=1;
    nm_per_pixel=1;
   
    switch nargin
        case 4
            imsize=varargin{4};       
        case 5           
            imsize=varargin{4};
            upFactor=varargin{5};
        case 6
            imsize=varargin{4};
            upFactor=varargin{5};
            nm_per_pixel=varargin{6};
        case 7
            imsize=varargin{4};
            upFactor=varargin{5};
            nm_per_pixel=varargin{6};
            aScale=varargin{7};
    end
end
    
    



%Converts 3, 2D arrays of cooridnates into an 3D image psuedo binary image.
row=imsize(1);
col=imsize(2);
stack=imsize(3);
if ~any(size(x)==1)

  %Interpolate within points
     xc=[x',x(1,:)']';
     yc=[y',y(1,:)']';
     zc=[z',z(1,:)']';    
    xc=interp2(xc,3);
    yc=interp2(yc,3);
    zc=interp2(zc,3);
    xc(end,:)=[];
    yc(end,:)=[];
    zc(end,:)=[];
    if exist('aScale','var')
        aScale=[aScale',aScale(1,:)']';
        aScale=interp2(aScale,3);
        aScale(end,:)=[];
    else
        
    [~,aScale]=surfArea(xc,yc,zc);
    end
    
    if all(aScale==0)
        aScale=ones(size(aScale));
    end
    
aScale=reshape(aScale,numel(aScale),1);
aScale=aScale/max(aScale);

        verteciesM=[reshape(xc,numel(xc),1),reshape(yc,numel(yc),1),...
            reshape(zc,numel(zc),1)];
        verteciesM=verteciesM/nm_per_pixel;
else
    
    while all(membrane2==Mfill)
F=triinterp2(F);
    membrane2=coord2image3d(F...
        ,imsize,upFactor,1,1)>0;
    membrane2=imclose(membrane2,true(5,5,5);
    Mfill=imfill(membrane2,'holes');
    end
    
    Mfill=Mfill-membrane2;
    [mx,my,mz]=find(membrane2);
    dist=pdist2([mx,my,mz],F.vertices,'euclidian','smallest');
    Mfill(membrane2)=1-dist;
    
%         x=interp1(x,1:1/10:numel(z))';
%     y=interp1(y,1:1/10:numel(z))';
%     z=interp1(z,1:1/10:numel(z))';
    verteciesM=[x,y,z];
    verteciesM=verteciesM/nm_per_pixel;
 %   [~,aScale]=gradient(verteciesM);
  %  aScale=sqrt(sum(aScale.^2,2))/nm_per_pixel;
end

Xpts=(0:1/upFactor:row)+1/upFactor/2;
Ypts=(0:1/upFactor:col)+1/upFactor/2;
Zpts=(0:1/upFactor:stack)+1/upFactor/2;



        [Image,~,~,locs]=histcn(verteciesM,Xpts,Ypts,Zpts);

    if (~any(size(x)==1) || exist('F','var'))

        for s=1:length(locs)
            if all(locs(s,:))    
            Image(locs(s,1),locs(s,2),locs(s,3))=(Image(locs(s,1),...
                locs(s,2),locs(s,3))-(1-aScale(s)));
            end
        end
    end