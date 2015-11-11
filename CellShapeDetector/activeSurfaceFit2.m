
function [x,y,z,eg]=activeSurfaceFit(x1,y1,z1,hyper_stack,cline_para,image_para,upFactor,flag)
minIterations=min(11,round(cline_para.iterations3d/2));
row=image_para.imsize(1);
col=image_para.imsize(2);
stack_z_size=image_para.imsize(3);
imsize=image_para.imsize;
ZX_factor=image_para.ZX_factor;
options.FourierKernel=true;
options.Power2Flag=false;
noiseIm=1./sqrt(hyper_stack+.01);
noiseIm=noiseIm./mean(noiseIm(hyper_stack>.4));
noiseIm=1;
%% Load PSF
averageSurfPSF = getappdata(0,'averageSurfPSF');
upFactor=2;
if flag.psf4d==1
    aPSF2=normalizeRange(averageSurfPSF);
    aPSF=aPSF2.^(1/upFactor);
    for i=1:size(aPSF,4)
        aPSF(:,:,:,i)=image_resize(aPSF2(:,:,:,i),size(aPSF2,1),...
            size(aPSF2,2),round(ZX_factor*size(aPSF2,3)));
    end
else
    aPSF=normalizeRange(averageSurfPSF(:,:,:,1));
    aPSF=image_resize(aPSF,size(aPSF,1),...
        size(aPSF,2),round(ZX_factor*size(aPSF,3)));
end

aPSF2=aPSF;
aPSF=image_resize(aPSF,size(aPSF,1)*upFactor,size(aPSF,2)*upFactor,...
    size(aPSF,3)*upFactor);
%aPSF=padarray(aPSF,[1,1,1]);

% aPSF2=flipdim(aPSF2,1);
% aPSF2=flipdim(aPSF2,2);
% aPSF2=aPSF2/2+flipdim(aPSF2,3)/2;
aPSFfft=fftn(aPSF,size(aPSF)+imsize*upFactor-[1,1,1]);

[yPSF,xPSF,zPSF]=gradient(aPSF2,1);
aPSFfft2=fftn(aPSF2,size(aPSF2)+imsize-[1,1,1]);
yPSFfft=fftn(yPSF,size(aPSF2)+imsize-[1,1,1]);
xPSFfft=fftn(xPSF,size(aPSF2)+imsize-[1,1,1]);
zPSFfft=fftn(zPSF,size(aPSF2)+imsize-[1,1,1]);



%% Create initial membrane outline and simulated image


% end_slices=max(end_slices,start_slices);
% r_vector=cellShellVectors(x1,y1,z1,end_slices);
[~,vertAreaI,V]=surfArea(x1,y1,z1);
V1=V;
V1=reshape(V1,[],3);
  [Xu(:,:,1),Xv(:,:,1)]=gradient(x1);
    [Xu(:,:,2),Xv(:,:,2)]=gradient(y1);
    [Xu(:,:,3),Xv(:,:,3)]=gradient(z1);
Xv1=reshape(Xv,numel(x1),3);
Xu1=reshape(Xu,numel(x1),3);
Xv1=normr(Xv1);
Xu1=normr(Xu1);
membrane=coord2image3d(x1,y1,z1,imsize,upFactor,1,vertAreaI);


x=x1;
y=y1;
z=z1;
vertices=[reshape(x1,numel(x1),1),reshape(y1,numel(y1),1),...
    reshape(z1,numel(z1),1)];

if cline_para.gradflag
membrane=fillMembrane(membrane,x,y,z,upFactor);
end

%vertices=vertices;

if flag.fftconv==1
    which convnfft
        convFit=convnfft(membrane, aPSFfft,'same',1:3,options);


else
    convFit=conv3sparse(membrane, aPSF,'same');
end

%% 3D Forces
%Calculate

%FL and FU are the upper and lower triangular factorizations that can
%be inverted
[FL,FU]=Finternal3D(x1,y1,z1,cline_para,'distance');
%cline_para.alpha =cline_para.stiff3d;
if ~isfield(cline_para,'zalpha3d')
cline_para.zalpha3d=cline_para.alpha3d;
cline_para.zbeta3d =cline_para.beta3d;
end
[FLz,FUz]=Finternal3Dz(x1,y1,z1,cline_para,'distance');



%% Cut off points that lie outside the image, scale simulated image
% to be similar to actual image
Xpts=(1:upFactor:upFactor*row)+1;
Ypts=(1:upFactor:upFactor*col)+1;
Zpts=(1:upFactor:upFactor*stack_z_size)+1;
[Ypts,Xpts,Zpts]=meshgrid(Ypts,Xpts,Zpts);

Xpts2=(0:1/upFactor:row)+1/upFactor/2;
Ypts2=(0:1/upFactor:col)+1/upFactor/2;
Zpts2=(0:1/upFactor:stack_z_size)+1/upFactor/2;
[Ypts2,Xpts2,Zpts2]=meshgrid(Ypts2,Xpts2,Zpts2);

%%

convFitSub=interp3(convFit,Ypts,Xpts,Zpts);

vertices(vertices(:,1)>row,1)=row;
vertices(vertices(:,2)>col,2)=col;
vertices(vertices(:,3)>stack_z_size,3)=stack_z_size;
vertices(vertices(:,1)<1,1)=1;
vertices(vertices(:,2)<1,2)=1;
vertices(vertices(:,3)<1,3)=1;



% imAvg=mean(hyper_stack(sub2ind(size(hyper_stack),round(vertices(:,1)),...
%     round(vertices(:,2)),round(vertices(:,3)))));
% fitAvg=mean(convFitSub(sub2ind(size(hyper_stack),round(vertices(:,1)),...
%     round(vertices(:,2)),round(vertices(:,3)))));
%         imAvg=mean(hyper_stack(:));
%         fitAvg=mean(convFitSub(:));
%K=mean(convFitSub(M2).*hyper_stack(M2))/mean(convFitSub(M2).^2);
K=mean(convFitSub(:).*hyper_stack(:))/mean(convFitSub(:).^2);


Kr=1;


deltaMembrane=zeros(size(membrane));
%membrane=membrane*K;
convFitSub=convFitSub*K;
convFit=convFit*K;

%% Update fits
vertices_avg=zeros(size(vertices));
best=Inf;
counter=0;
for i=1:cline_para.iterations3d
    membrane=membrane+deltaMembrane;
    if flag.fftconv==1
        convFit=convnfft(membrane, aPSFfft,'same',1:3,options);      
    else
        deltaConv=conv3sparse(deltaMembrane,aPSF,'same');
        convFit=convFit+deltaConv;
    end
        convFitSub=interp3(convFit,Ypts,Xpts,Zpts);

    %scale simulated image to be similar to actual image
%     imAvg=mean(hyper_stack(sub2ind(size(hyper_stack),round(vertices(:,1)),...
%         round(vertices(:,2)),round(vertices(:,3)))));
%     fitAvg=mean(convFitSub(sub2ind(size(hyper_stack),round(vertices(:,1)),...
%         round(vertices(:,2)),round(vertices(:,3)))));
%     
    
    if 1
        %scale image to minimize energy
        Kr=1;

        %Kr2=mean(convFitSub(M2).*hyper_stack(M2))/mean(convFitSub(M2).^2);
        Kr2=mean(convFitSub(:).*hyper_stack(:))/mean(convFitSub(:).^2);

    else
       %scale pointes individually to make simulated point more similar
       %no longer working
    Kr2=imAvg/fitAvg;
    Kr=hyper_stack./convFitSub;
    Kr(convFitSub<.1)=1;
    Kr=interp3(Kr,Ypts2,Xpts2,Zpts2);
    Kr(isnan(Kr) | isinf(Kr))=1;
    end
   % K=K.*Kr;
   % membrane=membrane.*Kr;
    convFit=convFit.*Kr2;
    convFitSub=convFitSub*Kr2;
    deltaIm=(convFitSub-hyper_stack);  % fix z scaling
    
    if flag.fftconv==0
        if i==1
            deltaImE=conv3sparse(deltaIm.*(convFitSub~=0),aPSF2,'same');
            deltaImOld=zeros(size(deltaImE));
        end
        deltaImD=deltaIm-deltaImOld;
        deltaImE=deltaImE+ conv3sparse(deltaImD,aPSF2,'same');%./areaCorrection;
        deltaImOld=deltaIm;
        
    end

        
    if cline_para.gradflag

            [~,~,V]=surfArea(x,y,z);
        V=reshape(V,[],3);
       
            deltaImConv=convnfft(deltaIm,aPSFfft2,'same',1:3,options);
    
    fs=interp3(deltaImConv,vertices(:,2),vertices(:,1),vertices(:,3),'*linear');    
        fs=bsxfun(@times,V,fs)/100;
    else
    if i<100
%         deltaImE=deltaIm;
%         [Fy,Fx,Fz]=(gradient(deltaIm,3));
                [Fy,Fx,Fz]=(gradient(convnfft(deltaIm,aPSFfft2,'same',1:3,options)));
    else
       
        [Fy,Fx,Fz]=(gradient(convnfft(noiseIm.*deltaIm,aPSFfft2,'same',1:3,options)));
    end
        fs=[interp3(Fx,vertices(:,2),vertices(:,1),vertices(:,3),'*linear'),...
        interp3(Fy,vertices(:,2),vertices(:,1),vertices(:,3),'*linear'),...
        interp3(Fz,vertices(:,2),vertices(:,1),vertices(:,3),'*linear')];
    end

    if i==minIterations
        eg=convFitSub;
        bestx=x;
        besty=y;
        bestz=z;
        counter=0;
        [~,vertAreaI,V]=surfArea(x,y,z);

    end
     if i==(minIterations-1);
         
         [~,vertAreaI,V]=surfArea(x,y,z);
     end
    
    dK=interp3(deltaIm,log2(upFactor));

    if (sum(dK(:).^2)<best&& i>minIterations) ||i==minIterations
        eg=convFitSub;
        bestx=x;
        besty=y;
        bestz=z;
        counter=0;
        best=sum(dK(:).^2);
    elseif i>minIterations
        counter=counter+1;
    end
    if counter==20;
        break
    end
    
    
    
    %Find distance to move by matrix division
    new_vertices(:,1:2) = FU\(FL \...
        (cline_para.gamma*vertices(:,1:2)-cline_para.kappa3d*fs(:,1:2)));
    new_vertices(:,3) = FUz\(FLz \ ...
        (cline_para.gamma*vertices(:,3)-cline_para.kappa3d*fs(:,3)));
    
    %if something goes wrong, backtrack
    if  all(isnan(new_vertices(:)))||~all(new_vertices(:)>.5)||...
            ~all(new_vertices(:,1)<row)||~all(new_vertices(:,2)<col)||...
            ~all(new_vertices(:,3)<stack_z_size)
        new_vertices(new_vertices<.5)=1;
        new_vertices(:,1)=min(new_vertices(:,1),row);
        new_vertices(:,2)=min(new_vertices(:,2),col);
        new_vertices(:,3)=min(new_vertices(:,3),stack_z_size);
        

    end

    deltaVert=new_vertices-vertices;
    
    if i==4
                  [Xu(:,:,1),Xv(:,:,1)]=gradient(x);
    [Xu(:,:,2),Xv(:,:,2)]=gradient(y);
    [Xu(:,:,3),Xv(:,:,3)]=gradient(z);
Xv1=reshape(Xv,numel(x1),3);
Xu1=reshape(Xu,numel(x1),3);
Xv1=normr(Xv1);
Xu1=normr(Xu1);
    end
    
    
    if i>4


circVec=bsxfun(@times,Xv1,dot((deltaVert),Xv1,2))+bsxfun(@times,Xu1,dot((deltaVert),Xu1,2));
    deltaVert(cline_para.pts+1:end-cline_para.pts,:)=...
        deltaVert(cline_para.pts+1:end-cline_para.pts,:)-...
        circVec(cline_para.pts+1:end-cline_para.pts,:);
    end
    vertices=vertices+deltaVert;

    vertices(1:cline_para.pts,:)=repmat(mean(vertices(1:cline_para.pts,:),1),cline_para.pts,1);
    vertices(end-cline_para.pts+1:end,:)=repmat(mean(vertices(end-cline_para.pts+1:end,:),1),cline_para.pts,1);
    %  Average over the last 20 iterations
%     if i>cline_para.iterations3d-20+1;
%         vertices_avg=vertices+vertices_avg;
%     end
%     if i==cline_para.iterations3d
%         vertices_avg=vertices+vertices_avg;
%         
%         vertices=vertices_avg/20;
%     end
    vertices2=reshape(vertices,numel(x1),3);
    x=reshape(vertices2(:,1),size(x1));
    y=reshape(vertices2(:,2),size(x1));
    z=reshape(vertices2(:,3),size(x1));
% Use these if you want to plot the surface with the force vecctors.
%     fsMatrix=reshape(fs,size(x1,1),size(x1,2),3);
%     deltaMatrix=reshape(new_vertices-vertices,size(x1,1),...
%         size(x1,2),3);
%     vertices2=reshape(vertices,numel(x1),3);
    
    membrane2=coord2image3d(x,y,z,imsize,upFactor,1,vertAreaI);
if cline_para.gradflag
membrane2=fillMembrane(membrane2,x,y,z,upFactor);
end


    
    if cline_para.show_flag>0 
        surf(x,y,z,vertAreaI);axis equal
        hold on
        
        quiver3(vertices(:,1),vertices(:,2),vertices(:,3),deltaVert(:,1),deltaVert(:,2),deltaVert(:,3));
        pause(.05)
        hold off
        [sum(dK(:).^2),best,i]
    end
  %  membrane2=membrane2.*K;
    
    deltaMembrane=membrane2-membrane;
    
end

x=bestx;
y=besty;
z=bestz;
end


