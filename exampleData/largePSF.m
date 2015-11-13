load dvpsf20120730.mat;
xtrim=15;
ytrim=15;
ztrim=10;

   xtrim=10;
   ytrim=10;
ztrim=0;

averageSurfPSF=normalizeRange(averagedPSF(100).addedPSF);
PSFsum=sum(averageSurfPSF(:));
averageSurfPSF=averageSurfPSF(xtrim+1:end-xtrim,ytrim+1:end-ytrim,ztrim+1:end-ztrim);
averageSurfPSF=averageSurfPSF/sum(averageSurfPSF(:))*PSFsum;
 averageSurfPSF=averageSurfPSF.^1;
 
  averageSurfPSF=(averageSurfPSF+flipdim(averageSurfPSF,1))/2;
  averageSurfPSF=(averageSurfPSF+flipdim(averageSurfPSF,2))/2;
averageSurfPSF=(averageSurfPSF+flipdim(averageSurfPSF,3))/2;
% 
% 
%averageSurfPSF(averageSurfPSF<.05)=0;
averageSurfPSF=normalizeRange(averageSurfPSF);
% mu = [0 0 0];
% Sigma = eye(3)*5;
% x1 = -15:15; x2 = -15:15; x3=-15:15;
% [X1,X2,X3] = meshgrid(x1,x2,x3);
% F = mvnpdf([X1(:) X2(:) X3(:)],mu,Sigma);
%averageSurfPSF=reshape(F,size(X1));
 setappdata(0,'averageSurfPSF',averageSurfPSF);
%clear all