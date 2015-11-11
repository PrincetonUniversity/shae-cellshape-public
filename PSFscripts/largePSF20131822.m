
load blockPSF20131822.mat

   xtrim=10;
   ytrim=10;
   ztrim=20;
aPSF=(averagedPSF(17).addedPSF);
aPSF(isnan(aPSF))=0;


%PSFsum=nansum(averageSurfPSF(:));
aPSF=aPSF(xtrim+1:end-xtrim,ytrim+1:end-ytrim,ztrim+1:end-ztrim);
aPSF(aPSF<.0)=0;
averageSurfPSF=normalizeRange(aPSF);
% mu = [0 0 0];
% Sigma = eye(3)*5;
% x1 = -15:15; x2 = -15:15; x3=-15:15;
% [X1,X2,X3] = meshgrid(x1,x2,x3);
% F = mvnpdf([X1(:) X2(:) X3(:)],mu,Sigma);
%averageSurfPSF=reshape(F,size(X1));
 setappdata(0,'averageSurfPSF',averageSurfPSF);
%clear all