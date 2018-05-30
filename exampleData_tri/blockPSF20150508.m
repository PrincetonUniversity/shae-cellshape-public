%%
load('averagedPSF20150508.mat')
%%
%    xtrim=10;
%    ytrim=10;
%    ztrim=20;

psf20150508 = averagedPSF(10).addedPSF./averagedPSF(10).sumMask;
psf20150508(isnan(psf20150508))=0;
psf20150508=psf20150508./max(psf20150508(:));
aPSF=psf20150508;
%PSFsum=nansum(averageSurfPSF(:));
% aPSF=aPSF(xtrim+1:end-xtrim,ytrim+1:end-ytrim,ztrim+1:end-ztrim);
%averageSurfPSF=avera[geSurfPSF/nansum(averageSurfPSF(:))*PSFsum;
%   aPSF=(aPSF+flipdim(aPSF,1))/2;
%   aPSF=(aPSF+flipdim(aPSF,2))/2;
% aPSF=(aPSF+permute(aPSF,[2,1,3]))/2;
  % averageSurfPSF=(averageSurfPSF+flipdim(averageSurfPSF,3))/2;
%  
% 
%averageSurfPSF(averageSurfPSF<.05)=0;
averageSurfPSF=normalizeRange(aPSF);
% mu = [0 0 0];
% Sigma = eye(3)*5;
% x1 = -15:15; x2 = -15:15; x3=-15:15;
% [X1,X2,X3] = meshgrid(x1,x2,x3);
% F = mvnpdf([X1(:) X2(:) X3(:)],mu,Sigma);
%averageSurfPSF=reshape(F,size(X1));
 setappdata(0,'averageSurfPSF',averageSurfPSF);
%clear all