
function [x,y,z,eg]=activeSurfaceFit_v1(x1,y1,z1,hyper_stack,cline_para,image_para,upFactor,flag)
row=image_para.imsize(1);
col=image_para.imsize(2);
stack_z_size=image_para.imsize(3);
imsize=image_para.imsize;
ZX_factor=image_para.ZX_factor;
options.FourierKernel=true;
options.Power2Flag=false;
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
aPSF=interpn(aPSF,log2(upFactor));
%aPSF=padarray(aPSF,[1,1,1]);

aPSF2=flipdim(aPSF2,1);
aPSF2=flipdim(aPSF2,2);
aPSF2=aPSF2/2+flipdim(aPSF2,3)/2;
aPSFfft=fftn(aPSF,size(aPSF)+imsize*upFactor-[1,1,1]);
[yPSF,xPSF,zPSF]=gradient(aPSF2,3);
yPSFfft=fftn(yPSF,size(aPSF2)+imsize-[1,1,1]);

xPSFfft=fftn(xPSF,size(aPSF2)+imsize-[1,1,1]);
zPSFfft=fftn(zPSF,size(aPSF2)+imsize-[1,1,1]);



%% Create initial membrane outline and simulated image


% end_slices=max(end_slices,start_slices);
% r_vector=cellShellVectors(x1,y1,z1,end_slices);

membrane=coord2image3d(x1,y1,z1,imsize,upFactor,1);

if flag.gradient
    membrane=imfill(double(membrane~=0),'holes');
end


vertecies=[reshape(x1,numel(x1),1),reshape(y1,numel(y1),1),...
    reshape(z1,numel(z1),1)];
%vertecies=vertecies;

if flag.fftconv==1
    which convnfft
        convFit=convnfft(membrane, aPSFfft,'same',1:3,options);


else
    convFit=conv3sparse(membrane, aPSF,'same');
end

%% 3D Forces
%Calculate
cline_para.alpha =cline_para.stiff3d;
cline_para.beta =cline_para.pts*cline_para.stiff3d;
%FL and FU are the upper and lower triangular factorizations that can
%be inverted
[FL,FU]=Finternal3D(x1,cline_para);
cline_para.alpha =cline_para.stiff3d;
cline_para.beta =4*cline_para.pts*cline_para.stiff3d;
[FLz,FUz]=Finternal3Dz(x1,cline_para);



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

vertecies(vertecies(:,1)>row,1)=row;
vertecies(vertecies(:,2)>col,2)=col;
vertecies(vertecies(:,3)>stack_z_size,3)=stack_z_size;
vertecies(vertecies(:,1)<1,1)=1;
vertecies(vertecies(:,2)<1,2)=1;
vertecies(vertecies(:,3)<1,3)=1;



imAvg=mean(hyper_stack(sub2ind(size(hyper_stack),round(vertecies(:,1)),...
    round(vertecies(:,2)),round(vertecies(:,3)))));
fitAvg=mean(convFitSub(sub2ind(size(hyper_stack),round(vertecies(:,1)),...
    round(vertecies(:,2)),round(vertecies(:,3)))));
K=imAvg/fitAvg;
if isnan(K)
    K=1;
end

Kr=1;


deltaMembrane=zeros(size(membrane));
membrane=membrane*K;
convFitSub=convFitSub*K;
convFit=convFit*K;
%% Update fits
vertecies_avg=zeros(size(vertecies));
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
    %scale simulated image to be similar to actual image
    imAvg=mean(hyper_stack(sub2ind(size(hyper_stack),round(vertecies(:,1)),...
        round(vertecies(:,2)),round(vertecies(:,3)))));
    fitAvg=mean(convFitSub(sub2ind(size(hyper_stack),round(vertecies(:,1)),...
        round(vertecies(:,2)),round(vertecies(:,3)))));
    if i>300
        Kr=1;
    else
        
    Kr2=imAvg/fitAvg;
    Kr=hyper_stack./convFitSub;
    Kr(convFitSub<.1)=1;
    Kr=interp3(Kr,Ypts2,Xpts2,Zpts2);
    Kr(isnan(Kr) | isinf(Kr))=1;
    end
    K=K.*Kr;
    membrane=membrane.*Kr;
    convFit=convFit.*Kr2;
    convFitSub=interp3(convFit,Ypts,Xpts,Zpts);
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
    
    
    if i<0
        deltaImE=deltaIm;
        %[Fy,Fx,Fz]=gradient(convnfft(deltaIm/100,aPSFfft,'same',1:3,options),3);
        [Fy,Fx,Fz]=(gradient(deltaIm,3));
    else
        
        Fz=convnfft(deltaIm,zPSFfft,'same',1:3,options)/100;
        Fx=convnfft(deltaIm,xPSFfft,'same',1:3,options)/100;
        Fy=convnfft(deltaIm,yPSFfft,'same',1:3,options)/100;
        
    end

    if i==11
        eg=convFitSub;
        bestx=x;
        besty=y;
        bestz=z;
        counter=0;
    end
    
    dK=interp3(deltaIm,log2(upFactor));
            [sum(dK(:).^2),Kr2,i]

    if (sum(dK(:).^2)<best&& i>11) ||i==11
        eg=convFitSub;
        bestx=x;
        besty=y;
        bestz=z;
        counter=0;
        best=sum(dK(:).^2);
    elseif i>11
        counter=counter+1;
    end
    if counter==15;
        break
    end
    
    
    %calculate forces at each point by interpolating
    fs=[interp3(Fx,vertecies(:,2),vertecies(:,1),vertecies(:,3),'*linear'),...
        interp3(Fy,vertecies(:,2),vertecies(:,1),vertecies(:,3),'*linear'),...
        interp3(Fz,vertecies(:,2),vertecies(:,1),vertecies(:,3),'*linear')];
    
    %Find distance to move by matrix division
    new_vertecies(:,1:2) = FU\(FL \...
        (cline_para.gamma*vertecies(:,1:2)-5*cline_para.kappa3d*fs(:,1:2)));
    new_vertecies(:,3) = FUz\(FLz \ ...
        (cline_para.gamma*vertecies(:,3)-5*cline_para.kappa3d*fs(:,3)));
    
    %if something goes wrong, backtrack
    if  all(isnan(new_vertecies(:)))||~all(new_vertecies(:)>.5)||...
            ~all(new_vertecies(:,1)<row)||~all(new_vertecies(:,2)<col)||...
            ~all(new_vertecies(:,3)<stack_z_size)
        new_vertecies(new_vertecies<.5)=1;
        new_vertecies(:,1)=min(new_vertecies(:,1),row);
        new_vertecies(:,2)=min(new_vertecies(:,2),col);
        new_vertecies(:,3)=min(new_vertecies(:,3),stack_z_size);
        

    end
    %After a few iterations, confine the points to move radially
    if i<8
        vertecies=new_vertecies;
        if i>1
%        r_vector=cellShellVectors(x,y,z,end_slices);
        end
    else
        vertecies=new_vertecies;%vertecies+(new_vertecies-vertecies).*r_vector;
    end
    
    vertecies(1:40,:)=repmat(mean(vertecies(1:40,:),1),40,1);
    vertecies(end-40+1:end,:)=repmat(mean(vertecies(end-40+1:end,:),1),40,1);
    %  Average over the last 20 iterations
    if i>cline_para.iterations3d-20+1;
        vertecies_avg=vertecies+vertecies_avg;
    end
    if i==cline_para.iterations3d
        vertecies_avg=vertecies+vertecies_avg;
        
        vertecies=vertecies_avg/20;
    end
    
    vertecies2=reshape(vertecies,numel(x1),3);
    x=reshape(vertecies2(:,1),size(x1));
    y=reshape(vertecies2(:,2),size(x1));
    z=reshape(vertecies2(:,3),size(x1));
% Use these if you want to plot the surface with the force vecctors.
%     fsMatrix=reshape(fs,size(x1,1),size(x1,2),3);
%     deltaMatrix=reshape(new_vertecies-vertecies,size(x1,1),...
%         size(x1,2),3);
%     vertecies2=reshape(vertecies,numel(x1),3);
    
    membrane2=coord2image3d(x,y,z,imsize,upFactor,1);
if flag.gradient
    membrane2=imfill(double(membrane2~=0),'holes');
end


    
    if cline_para.show_flag>0 
        surf(x,y,z);axis equal
        pause(.05)
    end

    membrane2=membrane2.*K;
    
    deltaMembrane=membrane2-membrane;
    sum(sum(sum(abs(deltaMembrane))));
    
end

x=bestx;
y=besty;
z=bestz;



% [kx,ky,kz]=meshgrid(1:201,1:201,1:201);
% kx=kx-101;
% ky=ky-101;
% kz=kz-101;
%K=(kx.^2+ky.^2+kz.^2);
%E=convnfft(membrane2,K,'same');
%E=normalizeRange(E);
%E=E.*mfill;
%E=E-min(E(E~=0));
%E=E+~mfill*1;
%E=E*100;






