%% RELATED TO DVPSF
% fit 4D DVPSF
% 
%%
averageSurfPSF = averagedPSF(100).addedPSF;
averageSurfPSF = averageSurfPSF./max(averageSurfPSF(:));

%% set up 'x' for 4D fitting
% list of xy pixels
X = 1:size(averageSurfPSF,2);
Y = 1:size(averageSurfPSF,1);
[xx2,yy2] = meshgrid(X-mean(X),Y-mean(Y));
% list of focal planes
zIDList = 10:22;

% list of positions from surface
zPrimeID = 80:100;
aOutStored = [];
% for zPrimeID = 90:100
zPSFsectionCenter = -5:0.05:1;
zPrimeList = zPSFsectionCenter(zPrimeID);

% starting point
a0 = [0.05,1,2,0,0,.1,0,1.5,0,0,0,0,0,0,0,0];
lb = ones(1,16).*-inf;
ub = ones(1,16).*inf;

% minimizeFunc
func2min = @(a) multiplaneComp(a,[xx2(:),yy2(:)],zIDList,zPrimeList,averagedPSF(zPrimeID));
options = optimset('MaxFun',3000,'MaxIter',500,'TolFun',1e-12);
[aOut,resnorm,residuals] = lsqnonlin(func2min,a0,lb,ub,options);

errorPSF = residuals;
biggestErrApprox = sqrt(max(errorPSF(:)).^2+min(errorPSF(:)).^2);
avgErr = sqrt(mean(errorPSF(:).^2));
minErr = min(errorPSF(:));
maxErr = max(errorPSF(:));
disp([zPrimeID(1), zPrimeID(end), biggestErrApprox, avgErr, minErr, maxErr]);
aOutStored = cat(1,aOutStored,aOut);

%%
% 
% a0 = [aaPoly1/365,aaPoly2,0.05,aaPoly4,aaPoly5,aaPoly6,16];
% residuals = func2min(a0);
% min(residuals)
% max(residuals)
% %%
% SliceBrowser(reshape(residuals,[51,51,13]));