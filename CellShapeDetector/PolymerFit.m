function [ poly ] = PolymerFit(F1, xf,yf,zf,Ff,cline_paraP,image_para )
%UNTITLED2 takes fluorescent stack image and cell coordinates x y and z,
%along with a fluorescent intensity on the surface, fitting parameters and
%image parameters
cline_paraP.alpha=0;
cline_paraP.beta=50000;
cline_paraP.gamma=40;
cline_paraP.kappa=100;
cline_paraP.stretching_force_factor=10;
cline_paraP.iterations=150;
cline_paraP.stretch_ends_flag=1;
cline_paraP.upFactor=1;
cline_paraP.thermal=0;
cline_paraP.objsize=3;
if ~isfield(cline_paraP,'hessFlag');
    cline_paraP.hessFlag=1;
end

%Ff2=interp3(bpass3_jn(F1,.25,[7,7,7]),xf,yf,zf);
%find spots on the unwrap, output is 3x wrapped
FZ=spotFind4(Ff,cline_paraP.hessFlag);
%clear objects near poles
FZ=padarray(FZ,[1,0],0);
FZ=FZ(:,2:end-1);
FZ=imclearborder(FZ);
FZ=FZ(2:end-1,:);
FZ=padarray(FZ,[0,1],0);

FZ=bwlabel(FZ>0,4);
clear rc

%3x repeat matrices
x3=[xf;xf;xf];
y3=[yf;yf;yf];
z3=[zf;zf;zf];
F3=[Ff;Ff;Ff];

meanI=accumarray(FZ(FZ>0),F3(FZ>0),[],@mean);
[meanI,ia]=sort(meanI,'descend');

FZ2=FZ;
for iblob=1:length(ia);
    FZ(FZ2==ia(iblob))=iblob;
end



%% make skeleton for each unique patch
for iblob=1:max(FZ(:))
FZmask=FZ==iblob;
FZmask=bwmorph(FZmask,'majority');
if all(any(FZmask'))
    FZmask([1:1:size(Ff,1),(2*size(Ff,1)+1):1:end],:)=0;
end

FZi=bwmorph(FZmask,'thin',Inf);
FZi=CleanSkeleton(FZi);
xyi=SkelOrder(FZi);
Yinit=xyi(:,1);
Xinit=xyi(:,2);
%[Yinit,Xinit]=find(FZi);

if numel(Xinit)==1
    Xinit=[Xinit(1);Xinit(end)+1];
    Yinit=[Yinit(1);Yinit(1)];
end
Xinit=smooth(Xinit,5);
Yinit=smooth(Yinit,5);

if length(Xinit)<10;
 Yinit=interp1(1/length(Xinit):1/length(Xinit):1,Yinit,1/10:1/10:1,'linear','extrap');
 Xinit=interp1(1/length(Xinit):1/length(Xinit):1,Xinit,1/10:1/10:1,'linear','extrap');
 cline_initial=[Xinit',Yinit'];

else
    cline_initial=[Xinit,Yinit];

end
cline_initial(:,~all(cline_initial))=[];

rc(iblob).coord=cline_initial;
end

%%

[poly,pts3d] = ActiveContourFitConvMembraneMulti3(F1, cline_paraP, rc,x3,y3,z3,F3,FZ);

for iPolymer=1:length(poly)
    poly(iPolymer).xs=poly(iPolymer).xs*image_para.nm_per_pixel;
    poly(iPolymer).ys=poly(iPolymer).ys*image_para.nm_per_pixel;
    poly(iPolymer).zs=poly(iPolymer).zs*image_para.nm_per_pixel;
    poly(iPolymer).rc=rc(iPolymer).coord;
    poly(iPolymer).length=poly(iPolymer).length*image_para.nm_per_pixel;
    poly(iPolymer).Curvature=poly(iPolymer).Curvature/image_para.nm_per_pixel;
    poly(iPolymer).Twist=poly(iPolymer).Twist*image_para.nm_per_pixel;
end

end

