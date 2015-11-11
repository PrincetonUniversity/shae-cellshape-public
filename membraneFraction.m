function [mFraction,Fmembrane] =membraneFraction(varargin)
% membraneFraction takes cell membrane coorrdinates and fluorescent image
% F1 with image parameters to estimate the fraction of fluorescent signal
% that is bound to membrane and the fraction that is cytoplasmic. Takes
% either patch object F or X Y and Z coordinates. PSF must be saved to the
% command window.

%usage
% mFraction=membraneFraction(tri,F1,image_para)
% or
% mFraction=membraneFraction(tri,F1)
%
%
% mFraction=membraneFraction(X,Y,Z,F1,image_para)
% or
% mFraction=membraneFraction(X,Y,Z,F1,image_para,hyper_stack,fillFlag)
% or
% mFraction=membraneFraction(X,Y,Z,F1)
%



if isfield(varargin{1},'faces')
    tri=varargin{1};
    F1=varargin{2};
    
    if nargin==2
        image_para.nm_per_pixel=1;
        image_para.imsize=size(F1);
        image_para.ZX_factor=.75;
    else
        image_para=varargin{3};
        
    end
    
    
    averageSurfPSF = getappdata(0,'averageSurfPSF');
    aPSF=normalizeRange(averageSurfPSF(:,:,:,1));
    aPSF=image_resize(aPSF,size(aPSF,1),...
        size(aPSF,2),round(image_para.ZX_factor*size(aPSF,3)));
    
    tri=triInterp2(tri,[],2);
    membrane2=coord2image3d(tri,image_para.imsize,1,image_para.nm_per_pixel);
    mfill=imfill(double(membrane2~=0),'holes');
    convMem=convnfft(membrane2,aPSF,'same');
    convMem=normalizeRange(convMem);
    convFill=convnfft(mfill,aPSF,'same');
    convFill=normalizeRange(convFill);
    
    x=fminsearch(@(x) sum(sum(sum((F1-(x(1)*convFill+(x(2))*convMem)).^2))),[.5,.5]);
    convFill=convFill.*x(1);
    Fmembrane=F1-convFill;
    Fmembrane(Fmembrane<0)=0;
    Ffill=max(sum(convFill(:)),0);
    Fmem=max(sum(Fmembrane(:)),0);
    mFraction=Fmem/(Fmem+Ffill);
    
else
    
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    F1=varargin{4};
    
    if nargin==5;%isempty(image_para)
        image_para.nm_per_pixel=1;
        image_para.imsize=size(F1);
        image_para.ZX_factor=.75;
    end
    if nargin>5
        image_para=varargin{5};
        hyper_stack=varargin{6};
        fillFlag=varargin{7};
    end
    
    
    averageSurfPSF = getappdata(0,'averageSurfPSF');
    aPSF=normalizeRange(averageSurfPSF(:,:,:,1));
    aPSF=image_resize(aPSF,size(aPSF,1),...
        size(aPSF,2),round(image_para.ZX_factor*size(aPSF,3)));
    
    
    membrane2=coord2image3d(X,Y,Z,image_para.imsize,1,image_para.nm_per_pixel);
    if ~nnz(membrane2(:))
        membrane2=coord2image3d(X,Y,Z,image_para.imsize,1,1);
    end
    
    mfill=imfill(double(membrane2~=0),'holes');
    
    if nargin>5
        if fillFlag
            convFill=normalizeRange(hyper_stack);
            convMem=convnfft(membrane2,aPSF,'same');
            convMem=normalizeRange(convMem);
        else
            convMem=normalizeRange(hyper_stack);
            convFill=convnfft(mfill,aPSF,'same');
            convFill=normalizeRange(convFill);
        end
    else
        convFill=convnfft(mfill,aPSF,'same');
        convFill=normalizeRange(convFill);
        convMem=convnfft(membrane2,aPSF,'same');
        convMem=normalizeRange(convMem);
    end
    
    
    [x,res]=fminsearch(@(x) sum(sum(sum((F1-(x(1)*convFill+(x(2))*convMem)).^2))),[.5,.5]);
    if x(2)<0
        [x,res]=fminsearch(@(x) sum(sum(sum((F1-(x(1)*convFill)).^2))),.5);
    end
    
    
    convFill=convFill.*x(1);
    Fmembrane=F1-convFill;
    Fmembrane(Fmembrane<0)=0;
    convFill=F1-Fmembrane;
    convFill(convFill<0)=0;
    Ffill=max(sum(convFill(:)),0);
    Fmem=max(sum(Fmembrane(:)),0);
    mFraction=Fmem/(Fmem+Ffill);
end