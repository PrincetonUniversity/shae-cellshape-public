%% group multiple files together 

% select a few files
whichFiles = uigetfile('*.mat','MultiSelect','on');
if not(iscell(whichFiles))
    whichFiles = {whichFiles};
end
savedPeaks2 = [];
for iFile = 1:size(whichFiles,2);
   load(whichFiles{iFile},'savedPeaks','firstStagePos','emFilter');
   switch emFilter
       case '534/30 nm (2)'
           savedPeaks(:,18) = 2;
       case '641/75 nm (3)'
           savedPeaks(:,18) = 3;
       otherwise
           savedPeaks(:,18) = nan;
   end
   savedPeaks(:,1) = savedPeaks(:,1)/1000+firstStagePos(1);
   savedPeaks(:,2) = savedPeaks(:,2)/1000+firstStagePos(2);
   savedPeaks(:,17) = savedPeaks(:,17)*(10^iFile);
   savedPeaks2 = cat(1,savedPeaks2,savedPeaks);
end

%%
hist(savedPeaks(:,3),100);
savedPeaks2(savedPeaks2(:,3)==0,:)=[];
newVar = [savedPeaks2(:,3),log(savedPeaks2(:,17))/6];
% newVar(newVar(:,2)<1.5,:) = [];

%% remove asymetric peaks
keepIdx = (savedPeaks2(:,14)./savedPeaks2(:,15))<1.15;
savedPeaks3 = savedPeaks2(keepIdx,:);
newVar = [savedPeaks3(:,3),log(savedPeaks3(:,17))/6];

[IDX,C] = kmeans(newVar,7,'replicates',50);

clf
hold on;
colrs = jet(max(IDX));
for ii = 1:max(IDX);
% for ii = 3;
scatter(newVar(IDX==ii,1),(newVar(IDX==ii,2)),15,colrs(ii,:),'o');
scatter(C(ii,1),C(ii,2),35,colrs(ii,:),'x','LineWidth',2);
end
axis equal;
figure(gcf);

%%
figure(gcf);
clf;
idx = savedPeaks2(:,18);
scatter3(savedPeaks2(idx==2,1),savedPeaks2(idx==2,2),savedPeaks2(idx==2,3),'gx');
hold on;
scatter3(savedPeaks2(idx==3,1),savedPeaks2(idx==3,2),savedPeaks2(idx==3,3),'rx');

%%

planeFunc = @(a,x) a(1)+a(2).*x(:,1)+a(3).*x(:,2);
    % fit to plane
    [bestPlaneOut,resnorm,residual] = lsqcurvefit(planeFunc,[1,1,1],...
        savedPeaks2(IDX==ii,1:2),savedPeaks2(IDX==ii,3));

    new3 = savedPeaks2;
    new3(:,3) = new3(:,3)-planeFunc(bestPlaneOut,new3(:,[1,2]));
    
    scatter(new3(:,3),log(new3(:,17)));