function FCorrect=FintensityAngleCorrection(x,y,z,F)
xc=bsxfun(@minus, x,mean(x));
yc=bsxfun(@minus, y,mean(y));
zc=bsxfun(@minus, z,mean(z));
rc=sqrt(xc.^2+yc.^2);
thetas=atan2(zc,rc);
thetas(thetas<-pi/2)=-pi/2+.01;
thetas(thetas>pi/2)=pi/2-.01;



[thetaCorrection,thetaLookup]=angleIntensityCorrection;
thetaCorrection=(thetaCorrection);
FCorrect=interp1(thetaLookup,thetaCorrection,thetas).*F;










