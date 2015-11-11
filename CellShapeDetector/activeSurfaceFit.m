
function [x,y,z,eg]=activeSurfaceFit(x1,y1,z1,hyper_stack,cline_para,image_para,upFactor,flag)
minIterations=min(11,round(cline_para.iterations3d/2));
VOXELS_COUNT_A = 15;

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
if isfield(cline_para,'upFactor');
    upFactor=cline_para.upFactor;
else
    upFactor=2;
end
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

%%%%%%%%%%%%%%
aPSF=pedistalSubtract(aPSF,'z2');
aPSF(aPSF<.0)=0;



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
% yPSFfft=fftn(yPSF,size(aPSF2)+imsize-[1,1,1]);
% xPSFfft=fftn(xPSF,size(aPSF2)+imsize-[1,1,1]);
% zPSFfft=fftn(zPSF,size(aPSF2)+imsize-[1,1,1]);



%% Create initial membrane outline and simulated image


% end_slices=max(end_slices,start_slices);
% r_vector=cellShellVectors(x1,y1,z1,end_slices);
[~,vertAreaI,V]=surfArea(x1,y1,z1);
V1=V;
V1=reshape(V1,[],3);
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

if flag.profile
    ticA = tic;
end

for i=1:cline_para.iterations3d
    
    if flag.profile
    if mod(i,10)==0;
        disp(['Time in active surface at i=',num2str(i),' is ', toc(ticA)]);
        toc(ticA);
    end
    end
        
    membrane=membrane+deltaMembrane;
    convFit=convnfft(membrane, aPSFfft,'same',1:3,options);
    convFitSub=interp3(convFit,Ypts,Xpts,Zpts);
    
    %scale simulated image to be similar to actual image
    %     imAvg=mean(hyper_stack(sub2ind(size(hyper_stack),round(vertices(:,1)),...
    %         round(vertices(:,2)),round(vertices(:,3)))));
    %     fitAvg=mean(convFitSub(sub2ind(size(hyper_stack),round(vertices(:,1)),...
    %         round(vertices(:,2)),round(vertices(:,3)))));
    %
    
    %scale image to minimize energy
    Kr=1;
    
    if cline_para.gradflag
        maxMask=hyper_stack>.8 & convFitSub>.5;
        approxObjSize=max(max(max(bwdist(~membrane))));
        if nnz(maxMask)<VOXELS_COUNT_A
            warning('activeSurfaceFit:voxelCountTooLow','Too few voxels above threshold for normalization');
            maxMask = hyper_stack>0.6 & convFitSub>0.4;
            if nnz(maxMask)<VOXELS_COUNT_A
                warning('activeSurfaceFit:voxelCountTooLow2','Too few voxels above secondary threshold for normalization');
                maxMask = hyper_stack>0.35 & convFitSub>0.2;
            end
        end
        % maxMask=imdilate(maxMask,true(3,3,3));
        maxMean=mean(convFitSub(maxMask));
        
        if approxObjSize>5
            convFitSub=(convFitSub.*(mean(hyper_stack(maxMask))))/maxMean*.95;
        else
            convFitSub=(convFitSub.*(mean(hyper_stack(maxMask))))/maxMean*.8;
        end
        
        
    else
        Kr2=mean(convFitSub(:).*hyper_stack(:))/mean(convFitSub(:).^2);
        convFitSub=convFitSub*Kr2;
        
    end
    
    % K=K.*Kr;
    % membrane=membrane.*Kr;
    deltaIm=(convFitSub-hyper_stack);  % fix z scaling
    
    
    if cline_para.gradflag
        
        [~,~,V]=surfArea(x,y,z);
        V=reshape(V,[],3);
        deltaIm=deltaIm.*~maxMask;%.*(hyper_stack>.2);
        
        deltaImConv=convnfft(deltaIm,aPSFfft2,'same',1:3,options);
        fs=interp3(deltaImConv,vertices(:,2),vertices(:,1),vertices(:,3),'*linear');
        fs2=reshape(fs,size(x));
        fs=bsxfun(@times,V,fs);
    else
        if i<100
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
    if  all(isnan(new_vertices(:)))||~all(new_vertices(:)>1)||...
            ~all(new_vertices(:,1)<row)||~all(new_vertices(:,2)<col)||...
            ~all(new_vertices(:,3)<stack_z_size)
        new_vertices(new_vertices<.5)=1;
        new_vertices(:,1)=min(new_vertices(:,1),row);
        new_vertices(:,2)=min(new_vertices(:,2),col);
        new_vertices(:,3)=min(new_vertices(:,3),stack_z_size);
        
        
    end
    
    deltaVert=new_vertices-vertices;
    if i>15
        deltaVert=bsxfun(@times,V1,dot((deltaVert),V1,2));
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
        if cline_para.gradflag
            surf(x,y,z,fs2);axis equal
        else
            surf(x,y,z);axis equal
        end
        
        colorbar
        hold on
        quiver3(vertices(:,1),vertices(:,2),vertices(:,3),fs(:,1),fs(:,2),fs(:,3));
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


